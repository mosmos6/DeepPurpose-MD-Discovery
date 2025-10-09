# 5_md_simulation.py  (MixMD-lite + JAX hotspots)
# Author: Iori Mochizuki + collab
# Updated: 2025-10-09
# Description: OpenMM-based MD (protein/RNA) with optional MixMD-lite.
# - --mixmd auto-enforces apo and adds small cosolvent probes (OpenFF)
# - JAX-based hotspot detection (grid occupancy, top-K peaks) from the produced DCD
# - Preserves existing filenames/flags so downstream scripts remain intact

import argparse, json, math, random, sys, os
from sys import stdout

# ---- OpenMM / OpenFF stack (unchanged core) ----
import openmm as mm
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.app import *
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---- Read-only trajectory load (MDTraj) ----
import mdtraj as md  # READ-ONLY: we only load xyz + cell; no MDTraj distance math

# ---- JAX: all distance/occupancy math done here (GPU-accelerated) ----
import jax
import jax.numpy as jnp
from jax import lax

# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand, no-ligand, MixMD-lite)")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Input ligand SDF")
parser.add_argument("--n-steps", type=int, default=1000000, help="Production MD steps (default: 1,000,000 ~2 ns)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# ---- NEW: MixMD-lite controls ----
parser.add_argument("--mixmd", action="store_true", help="Enable MixMD-lite (adds small probes; forces apo)")
parser.add_argument("--probe-total", type=int, default=150, help="Total number of cosolvent molecules to try to place")
parser.add_argument("--probe-fractions", type=str,
                    default="ipa:0.40,acn:0.20,imd:0.15,aceam:0.15,phol:0.10",
                    help="Comma-separated fractions (name:frac). Names: ipa,acn,imd,aceam,phol,acoh(optional)")
parser.add_argument("--grid-spacing", type=float, default=0.20, help="Hotspot grid spacing in nm (default 0.20 nm ~2 Å)")
parser.add_argument("--topk", type=int, default=3, help="Number of hotspot centers to report")
parser.add_argument("--analysis-skip-frames", type=int, default=0, help="Skip first N frames during occupancy (e.g., to skip early equilibration)")

args = parser.parse_args()

# MixMD implies apo
if args.mixmd:
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---------------------------
# Small helpers
# ---------------------------
def cosolvent_library():
    # neutral, chemically diverse set; acetate replaced by acetic acid (avoid charge juggling)
    # SMILES chosen for robust OpenFF parameterization
    return {
        "ipa":   {"name": "IPA",   "smiles": "CC(O)C"},          # isopropanol
        "acn":   {"name": "ACN",   "smiles": "CC#N"},            # acetonitrile
        "imd":   {"name": "IMD",   "smiles": "c1ncc[nH]1"},      # imidazole (neutral tautomer)
        "aceam": {"name": "ACEA",  "smiles": "CC(=O)N"},         # acetamide
        "phol":  {"name": "PHOL",  "smiles": "c1ccccc1O"},       # phenol
        "acoh":  {"name": "ACOH",  "smiles": "CC(=O)O"},         # acetic acid (optional)
    }

def parse_probe_fractions(s):
    lib = cosolvent_library()
    raw = [p.strip() for p in s.split(",") if p.strip()]
    kv = {}
    for piece in raw:
        k, v = piece.split(":")
        k = k.strip().lower()
        if k not in lib:
            raise ValueError(f"Unknown probe key '{k}'. Allowed: {','.join(lib.keys())}")
        kv[k] = float(v)
    total = sum(kv.values())
    if total <= 0:
        raise ValueError("Sum of probe fractions must be > 0.")
    # Normalize
    for k in kv:
        kv[k] = kv[k] / total
    return kv

def random_rotation_matrix(rng):
    # Uniform random rotation via quaternion method
    u1, u2, u3 = random.random(), random.random(), random.random()
    q1 = math.sqrt(1-u1) * math.sin(2*math.pi*u2)
    q2 = math.sqrt(1-u1) * math.cos(2*math.pi*u2)
    q3 = math.sqrt(u1)   * math.sin(2*math.pi*u3)
    q4 = math.sqrt(u1)   * math.cos(2*math.pi*u3)
    # Rotation matrix
    R = [[1-2*(q3*q3+q4*q4),   2*(q2*q3-q1*q4),     2*(q2*q4+q1*q3)],
         [2*(q2*q3+q1*q4),     1-2*(q2*q2+q4*q4),   2*(q3*q4-q1*q2)],
         [2*(q2*q4-q1*q3),     2*(q3*q4+q1*q2),     1-2*(q2*q2+q3*q3)]]
    return R

def to_jnp_positions(positions):
    """OpenMM Vec3 (nm) -> jnp array (N,3) in nm"""
    return jnp.array([[p.x, p.y, p.z] for p in positions], dtype=jnp.float32)

def minimum_image(delta, box_lengths):
    # delta, box_lengths are (N,3) and (3,) in nm
    return delta - jnp.round(delta / box_lengths) * box_lengths

def topk_peaks_from_grid(grid, topk, spacing, origin):
    # grid: jnp array (nx,ny,nz) counts
    flat = grid.reshape((-1,))
    vals, idx = lax.top_k(flat, topk)
    nx, ny, nz = grid.shape
    iz = idx % nz
    iy = (idx // nz) % ny
    ix = (idx // (ny * nz))
    centers_nm = jnp.stack([ix, iy, iz], axis=1).astype(jnp.float32) + 0.5
    centers_nm = origin + centers_nm * spacing
    return centers_nm, vals

def safe_int_tuple(x):
    return (int(x[0]), int(x[1]), int(x[2]))

# ---------------------------
# 1. Load receptor
# ---------------------------
print("📦 Loading receptor:", args.input_receptor)
receptor_pdb = PDBFile(args.input_receptor)

# ---------------------------
# 2. Force field setup
# ---------------------------
if args.rna:
    print("🧬 [RNA MODE] amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

# ---------------------------
# 3. Ligand logic (unchanged unless MixMD forces apo)
# ---------------------------
if not args.no_ligand and not args.mixmd:
    print("🔗 Holo mode with ligand:", args.input_ligand)
    ligand = Molecule.from_file(args.input_ligand)
    ligand_positions = to_openmm(ligand.conformers[0])
    ligand_top = OFFTopology.from_molecules([ligand]).to_openmm()
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[ligand])
    forcefield.registerTemplateGenerator(smirnoff.generator)
    modeller.add(ligand_top, ligand_positions)
else:
    print("🔗 Apo mode (no ligand).")

# ---------------------------
# 4. MixMD-lite probe placement (before solvation)
# ---------------------------
placed_probes = []
probe_defs = cosolvent_library()
probe_plan = {}

if args.mixmd:
    print("🧪 [MixMD-lite] Planning cosolvent placement...")
    fractions = parse_probe_fractions(args.probe_fractions)
    # compute counts per type
    for key, frac in fractions.items():
        probe_plan[key] = int(round(frac * args.probe_total))
    # Build OpenFF molecules & register for SMIRNOFF
    mols = []
    for key, cnt in probe_plan.items():
        if cnt <= 0: 
            continue
        info = probe_defs[key]
        m = Molecule.from_smiles(info["smiles"], allow_undefined_stereo=True)
        m.name = info["name"]
        m.generate_conformers(n_conformers=1)
        mols.append(m)
    if mols:
        smirnoff_mix = SMIRNOFFTemplateGenerator(molecules=mols)
        forcefield.registerTemplateGenerator(smirnoff_mix.generator)

    # Receptor centroid (for rough placement cube)
    rec_coords = to_jnp_positions(receptor_pdb.positions)
    rec_center = jnp.mean(rec_coords, axis=0)  # nm

    # Try to place probes uniformly in a ~6 nm cube around receptor center
    # (box to be defined by addSolvent later)
    rng = random.Random(args.seed)
    hard_min = 0.18  # nm min distance to receptor atoms (soft repel)
    half_box = 3.0   # nm half-side target for initial placement

    # Precompute receptor coords as plain lists for speed here
    rec_xyz = [[float(p.x), float(p.y), float(p.z)] for p in receptor_pdb.positions]

    total_target = sum(probe_plan.values())
    added = 0
    for key, cnt in probe_plan.items():
        if cnt <= 0: 
            continue
        info = probe_defs[key]
        # Prepare this probe's topology and conformer positions
        m = next(mm for mm in mols if mm.name == info["name"])
        top = OFFTopology.from_molecules([m]).to_openmm()
        base_pos = to_openmm(m.conformers[0])  # quantity Vec3 (Å internally -> OpenMM handles units)
        # Convert to numeric for transform (we’ll convert back)
        base_nm = jnp.array([[p.x, p.y, p.z] for p in base_pos], dtype=jnp.float32)  # should be nm already
        # Place cnt copies
        n_placed = 0
        trials = 0
        while n_placed < cnt and trials < cnt * 50:
            trials += 1
            # random rotation
            R = jnp.array(random_rotation_matrix(rng), dtype=jnp.float32)
            rot = (base_nm @ R.T)
            # random translation near receptor center
            t = jnp.array([
                float(rec_center[0]) + (rng.random()*2-1)*half_box,
                float(rec_center[1]) + (rng.random()*2-1)*half_box,
                float(rec_center[2]) + (rng.random()*2-1)*half_box,
            ], dtype=jnp.float32)
            cand = rot + t  # (Natoms,3) nm

            # quick hard-min check to receptor
            ok = True
            for i in range(cand.shape[0]):
                cx, cy, cz = float(cand[i,0]), float(cand[i,1]), float(cand[i,2])
                # sample a few receptor atoms (downsample for speed)
                for rx, ry, rz in rec_xyz[:: max(1, len(rec_xyz)//4000) ]:
                    dx = cx - rx; dy = cy - ry; dz = cz - rz
                    if (dx*dx + dy*dy + dz*dz) < (hard_min*hard_min):
                        ok = False; break
                if not ok:
                    break
            if not ok:
                continue

            # convert cand -> OpenMM Vec3 with nm units
            cand_pos = [mm.Vec3(float(cand[i,0]), float(cand[i,1]), float(cand[i,2])) for i in range(cand.shape[0])]
            modeller.add(top, cand_pos)
            placed_probes.append({"key": key, "name": info["name"]})
            n_placed += 1
            added += 1
        print(f"  • {info['name']}: placed {n_placed}/{cnt} (tried {trials})")

    # Save list for provenance
    with open("mixmd_probe_placement.json", "w") as f:
        json.dump({"plan": probe_plan, "placed": placed_probes}, f, indent=2)
    print(f"✅ MixMD probes placed: {added} total.")

# ---------------------------
# Hydrogens & write combined
# ---------------------------
modeller.addHydrogens(forcefield)
with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# ---------------------------
# Solvation (unchanged)
# ---------------------------
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(7.0, 7.0, 7.0) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("✅ Solvated system ready.")

# ---------------------------
# System creation (unchanged)
# ---------------------------
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds
)

seed = int(args.seed)
barostat = MonteCarloBarostat(1 * bar, 300 * kelvin, 25)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
integrator.setRandomNumberSeed(seed)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print("🔹 Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter(f"minimized{suffix}.pdb", 100))

print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * kelvin, seed)
simulation.reporters.append(PDBReporter(f"nvt_equilibrated{suffix}.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

print("🔹 NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter(f"npt_equilibrated{suffix}.pdb", 500))
simulation.step(2500)  # 5 ps

print(f"🔥 Production MD: {args.n_steps:,} steps")
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(args.n_steps)

state = simulation.context.getState(getPositions=True)
with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"🎉 MD complete → final_structure{suffix}.pdb saved.")

# ---------------------------
#  Hotspot analysis (JAX): MixMD only
# ---------------------------
if args.mixmd:
    print("🔎 [MixMD-lite] JAX hotspot analysis starting...")
    # Load trajectory (READ-ONLY via MDTraj)
    dcd = f"production_md{suffix}.dcd"
    top = f"solvated_receptor_ligand{suffix}.pdb"
    traj = md.load_dcd(dcd, top=top)  # nm units
    xyz = jnp.array(traj.xyz, dtype=jnp.float32)  # (F, A, 3) nm
    F, A, _ = xyz.shape

    # Identify selections using MDTraj topology (no distance calls)
    topol = traj.topology
    prot_idx = topol.select("protein")
    rna_idx = topol.select("nucleic")
    water_idx = topol.select("water")
    ion_idx = topol.select("name NA or name CL or name K or name MG or name CA")
    receptor_idx = jnp.array(sorted(set(prot_idx.tolist() + rna_idx.tolist())), dtype=jnp.int32)

    # Heavy atoms not receptor/water/ions -> "cosolvent pool"
    heavy_mask = jnp.array([0 if (a.element.symbol == 'H') else 1 for a in topol.atoms], dtype=jnp.int32)
    exclude = set(receptor_idx.tolist()) | set(water_idx.tolist()) | set(ion_idx.tolist())
    cosolv_list = [i for i in range(A) if (i not in exclude) and (heavy_mask[i] == 1)]
    cosolv_idx = jnp.array(cosolv_list, dtype=jnp.int32)

    if cosolv_idx.size == 0:
        print("⚠️ No cosolvent heavy atoms detected after solvation; skipping hotspot detection.")
        sys.exit(0)

    # Per-frame box vectors (nm)
    cell = jnp.array(traj.unitcell_lengths, dtype=jnp.float32)  # (F, 3) nm

    # Skip initial frames if requested
    f0 = max(0, int(args.analysis_skip_frames))
    xyz = xyz[f0:]
    cell = cell[f0:]
    F = xyz.shape[0]

    # Receptor COM per frame
    rec_xyz = xyz[:, receptor_idx, :]  # (F, Nr, 3)
    rec_com = jnp.mean(rec_xyz, axis=1)  # (F, 3)

    # Minimum-image unwrap of cosolvent coordinates around receptor COM
    cos_xyz = xyz[:, cosolv_idx, :]      # (F, Nc, 3)
    # Broadcasted minimum image
    # delta = r_cos - r_com
    delta = cos_xyz - rec_com[:, None, :]
    # Wrap per frame with that frame's box
    L = cell  # (F, 3)
    # Match shapes for broadcasting
    delta_mi = minimum_image(delta, L[:, None, :])
    cos_uw = rec_com[:, None, :] + delta_mi  # (F, Nc, 3)

    # Build a grid around receptor heavy atoms
    rec_min = jnp.min(rec_xyz.reshape(F, -1, 3), axis=(0,1))
    rec_max = jnp.max(rec_xyz.reshape(F, -1, 3), axis=(0,1))
    pad = jnp.array([0.5, 0.5, 0.5], dtype=jnp.float32)  # nm padding
    origin = rec_min - pad
    spacing = float(args.grid_spacing)
    # number of bins along each axis
    nxyz = jnp.ceil((rec_max + pad - origin) / spacing).astype(jnp.int32)
    nx, ny, nz = int(nxyz[0]), int(nxyz[1]), int(nxyz[2])
    print(f"🧱 Grid: {nx} x {ny} x {nz} (spacing={spacing:.2f} nm)")

    # Accumulate occupancy: flatten frames then bin
    # coords -> indices
    coords = cos_uw.reshape((-1, 3))  # (F*Nc, 3)
    rel = coords - origin[None, :]
    idxf = jnp.floor(rel / spacing).astype(jnp.int32)
    # clamp to grid
    idx0 = jnp.clip(idxf[:, 0], 0, nx-1)
    idx1 = jnp.clip(idxf[:, 1], 0, ny-1)
    idx2 = jnp.clip(idxf[:, 2], 0, nz-1)
    # flat indexing
    flat_idx = (idx0 * (ny * nz) + idx1 * nz + idx2).astype(jnp.int32)
    # scatter-add into flat grid
    grid_flat = jnp.zeros((nx * ny * nz,), dtype=jnp.int32)
    grid_flat = grid_flat.at[flat_idx].add(1)
    grid = grid_flat.reshape((nx, ny, nz))

    # Top-K peaks
    centers_nm, counts = topk_peaks_from_grid(grid, int(args.topk), spacing, origin)
    centers_A = centers_nm * 10.0  # nm -> Å

    # Write CSV in Å for Vina
    out_csv = "mixmd_hotspots.csv"
    with open(out_csv, "w") as f:
        f.write("rank,center_x,center_y,center_z,raw_count\n")
        for k in range(min(args.topk, centers_A.shape[0])):
            cx, cy, cz = float(centers_A[k,0]), float(centers_A[k,1]), float(centers_A[k,2])
            cnt = int(counts[k])
            f.write(f"{k+1},{cx:.3f},{cy:.3f},{cz:.3f},{cnt}\n")
    print(f"✅ Hotspots written -> {out_csv}")

    # Console helper for your Step 3
    if centers_A.shape[0] > 0:
        cx, cy, cz = float(centers_A[0,0]), float(centers_A[0,1]), float(centers_A[0,2])
        print("\n👉 Suggested Vina center (top hotspot):")
        print(f"   --center_x {cx:.3f} --center_y {cy:.3f} --center_z {cz:.3f}")
        print("   (Use with scripts/3_docking_vina.py; box size remains your default 30 Å.)")

print("✅ Done.")
