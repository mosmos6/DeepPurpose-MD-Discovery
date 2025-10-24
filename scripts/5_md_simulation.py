# 5_md_simulation.py
# Author: Iori Mochizuki (pipeline owner) + collaborator patch
# Updated: 2025-10-23
#
# OpenMM MD (protein/RNA/ligand) with:
#  - Multi-ligand support (auto-discovers ligand_*.sdf; falls back to ligand.sdf)
#  - MixMD-from-Packmol (adds cosolvent probes; logs probe centroids to CSV)
#  - Robust hotspot logging (segmented stepping; no reporter scheduling pitfalls)
#  - Auto mmCIF write for large systems (avoids PDB overflow warnings)
#
# File names remain as in your pipeline where possible:
#   combined_receptor_ligand[_no_ligand].pdb/.cif
#   solvated_receptor_ligand[_no_ligand].pdb/.cif
#   minimized[_no_ligand].pdb       nvt_equilibrated[_no_ligand].pdb
#   npt_equilibrated[_no_ligand].pdb production_md[_no_ligand].dcd/.log
#   final_structure[_no_ligand].pdb/.cif
#   mixmd_hotspots[_no_ligand].csv  (when --mixmd-from-packmol)
#
import argparse, csv, glob, os
from pathlib import Path
from sys import stdout

import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---------------------------
# MixMD probe library
# ---------------------------
PROBE_LIBRARY = {
    # canonical 3-letter → SMILES
    "IPA":  "CC(C)O",           # isopropanol
    "ACN":  "CC#N",             # acetonitrile
    "IMD":  "c1[nH]cnc1",       # imidazole (neutral, N1–H)
    "ACM":  "CC(=O)N",          # acetamide
    "PHO":  "c1ccc(cc1)O",      # phenol
    "HAC":  "CC(=O)O",          # acetic acid
    # legacy aliases (accepted on CLI)
    "ACEA": "CC(=O)N",
    "PHOL": "c1ccc(cc1)O",
    "ACOH": "CC(=O)O",
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())

# ---------------------------
# Small IO helpers
# ---------------------------
def _default_packmol_csv(receptor_path: Path) -> Path:
    root = f"{receptor_path.stem}_mixmd"
    return Path("build") / f"{root}_placements.csv"

def _write_structure_auto(topology, positions, stem: str) -> str:
    """
    Write PDB for smaller systems; mmCIF for very large (to avoid PDB overflow).
    Returns the written file path.
    """
    n_atoms = topology.getNumAtoms()
    n_res   = sum(1 for _ in topology.residues())
    # PDB column limits: 99999 atoms / 9999 residues
    if n_atoms > 99999 or n_res > 9999:
        out = f"{stem}.cif"
        with open(out, "w") as fh:
            PDBxFile.writeFile(topology, positions, fh)
    else:
        out = f"{stem}.pdb"
        with open(out, "w") as fh:
            PDBFile.writeFile(topology, positions, fh)
    return out

# ---------------------------
# Multi‑ligand discovery & merge
# ---------------------------
def _discover_ligand_sdfs() -> list[Path]:
    # prefer enumerated set ligand_*.sdf; fallback to single ligand.sdf
    paths = sorted(Path(".").glob("ligand_*.sdf"),
                   key=lambda p: (len(p.stem), p.stem))
    if paths:
        return paths
    if Path("ligand.sdf").exists():
        return [Path("ligand.sdf")]
    return []

def _merge_ligands_into_modeller(modeller: Modeller, forcefield: ForceField, sdf_paths: list[Path]) -> int:
    """
    Parameterize all ligands (SMIRNOFF) and append to modeller with their conformer positions.
    Returns number of ligands added.
    """
    if not sdf_paths:
        return 0
    off_mols = []
    for sdf in sdf_paths:
        m = Molecule.from_file(str(sdf))
        if not m.conformers:
            # ensure at least one conformer exists for placement
            m.generate_conformers(n_conformers=1)
        off_mols.append(m)
    # register one SMIRNOFF generator for all ligands present
    forcefield.registerTemplateGenerator(SMIRNOFFTemplateGenerator(molecules=off_mols).generator)

    added = 0
    for m in off_mols:
        top = OFFTopology.from_molecules([m]).to_openmm()
        pos = to_openmm(m.conformers[0])
        modeller.add(top, pos)
        added += 1
    return added

# ---------------------------
# Packmol → add probes
# ---------------------------
def _read_centroids_csv(csv_path: Path):
    rows = []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({
                "resname": row["resname"].strip().upper(),
                "resid":   int(row["resid"]),
                "x_nm":    float(row["x_nm"]),
                "y_nm":    float(row["y_nm"]),
                "z_nm":    float(row["z_nm"]),
            })
    return rows

def _openff_conf_to_angstrom_list(off_mol):
    conf = off_mol.conformers[0]
    try:
        arr = conf.to("angstrom").m  # Pint path
        return [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
    except Exception:
        q = to_openmm(conf)  # OpenMM Quantity(list(Vec3), angstrom)
        return [(float(v[0]), float(v[1]), float(v[2])) for v in q.value_in_unit(angstrom)]

def _positions_for_centroid_nm(off_mol, centroid_nm_tuple):
    coords_A = _openff_conf_to_angstrom_list(off_mol)
    n = len(coords_A)
    sx = sy = sz = 0.0
    for xA, yA, zA in coords_A:
        sx += xA; sy += yA; sz += zA
    cx = sx / n; cy = sy / n; cz = sz / n
    tx, ty, tz = float(centroid_nm_tuple[0]), float(centroid_nm_tuple[1]), float(centroid_nm_tuple[2])
    placed = []
    for xA, yA, zA in coords_A:
        xnm = (xA - cx) * 0.1 + tx
        ynm = (yA - cy) * 0.1 + ty
        znm = (zA - cz) * 0.1 + tz
        placed.append(mm.Vec3(xnm, ynm, znm))
    return [p * nanometer for p in placed]

def _add_probes_from_packmol(modeller: Modeller,
                             forcefield: ForceField,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None,
                             resname_list: list[str] | None) -> int:
    if placements_csv is None:
        placements_csv = _default_packmol_csv(receptor_path)
    print(f"ℹ️  MixMD inputs: CSV={placements_csv}")
    if not placements_csv.exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")

    half = 0.5 * float(sim_box_nm)
    keep_lo = - (half - float(edge_margin_nm))
    keep_hi = + (half - float(edge_margin_nm))
    whitelist = {r.upper().strip() for r in resname_list} if resname_list else None

    kept = []
    for rec in _read_centroids_csv(placements_csv):
        resname = rec["resname"]
        if resname not in PROBE_LIBRARY:
            continue
        if whitelist and resname not in whitelist:
            continue
        x, y, z = rec["x_nm"], rec["y_nm"], rec["z_nm"]
        if (keep_lo <= x <= keep_hi) and (keep_lo <= y <= keep_hi) and (keep_lo <= z <= keep_hi):
            kept.append((resname, rec["resid"], (x, y, z)))

    if not kept:
        print("⚠️  No probe centroids inside the MD box; continuing with receptor-only.")
        return 0

    used_res = sorted({res for (res, _, _) in kept})
    off_mols = []
    for res in used_res:
        m = Molecule.from_smiles(PROBE_LIBRARY[res], allow_undefined_stereo=True)
        m.generate_conformers(n_conformers=1)
        m.name = res
        off_mols.append(m)
    forcefield.registerTemplateGenerator(SMIRNOFFTemplateGenerator(molecules=off_mols).generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    # helper to see new residues for renaming
    def _res_list():
        return list(modeller.topology.residues())

    added = 0
    for (resname, pack_resid, cen_nm) in kept:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)

        before = len(_res_list())
        modeller.add(top, pos)
        tail = _res_list()[before:]

        for res in tail:
            try: res.name = resname
            except Exception: pass
            try: res.id = str(pack_resid)
            except Exception: pass

        added += 1

    print(f"✅ Placed {added} probe molecules (cropped to {sim_box_nm:.1f} nm box).")
    return added

# ---------------------------
# Hotspot logging (no reporter scheduling)
# ---------------------------
def _build_probe_groups(topology, allowed_names: set[str]):
    groups, resid_map = [], []
    for res in topology.residues():
        rn = (res.name or "").upper().strip()
        if rn in allowed_names:
            idxs = [a.index for a in res.atoms()]
            if idxs:
                groups.append(idxs)
                resid_map.append((rn, int(res.id) if res.id is not None else 0))
    return groups, resid_map

def _protein_anchor_indices(topology):
    aa = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
    }
    ca = []
    heavy = []
    for res in topology.residues():
        rn = (res.name or "").upper().strip()
        if rn in aa:
            for a in res.atoms():
                nm = (a.name or "").upper()
                if nm == "CA":
                    ca.append(a.index)
                if not (a.element and a.element.symbol == "H"):
                    heavy.append(a.index)
    return ca if ca else heavy

def _com_of(positions, idxs):
    sx = sy = sz = 0.0
    n = 0
    for i in idxs:
        v = positions[i]
        sx += float(v.x); sy += float(v.y); sz += float(v.z)
        n += 1
    if n == 0: return (0.0, 0.0, 0.0)
    return (sx/n, sy/n, sz/n)

def _append_hotspots_csv(fh, step, positions, groups, resid_map, recenter_dxdy dz):
    dx, dy, dz = recenter_dxdy
    for (idxs, (resname, resid)) in zip(groups, resid_map):
        sx = sy = sz = 0.0
        n = 0
        for i in idxs:
            v = positions[i]
            sx += float(v.x); sy += float(v.y); sz += float(v.z)
            n += 1
        if n > 0:
            fh.write(f"{step},{resname},{resid},{sx/n+dx:.5f},{sy/n+dy:.5f},{sz/n+dz:.5f}\n")

def _run_segmented(simulation: Simulation,
                   total_steps: int,
                   stride: int,
                   csv_path: Path | None,
                   probe_names: set[str] | None):
    """
    If csv_path is provided, log probe centroids every 'stride' steps by:
      - building groups once (from simulation.topology)
      - computing receptor COM(0), then at each checkpoint shifting coordinates by
        COM(0) - COM(t) (protein-centric frame).
    Returns the final absolute step reached (int).
    """
    step_counter = 0
    fh = None
    groups = resid_map = None
    anchor = None
    ref_com = None
    try:
        if csv_path:
            fh = open(csv_path, "w")
            names = probe_names or set()
            # autodiscover + user set
            present = {(r.name or "").upper().strip() for r in simulation.topology.residues()}
            allowed = (present & PROBE_NAME_SET) | (names & PROBE_NAME_SET)
            fh.write("# tracking_resnames=" + ";".join(sorted(allowed)) + "\n")
            fh.write("# center_mode=protein-com\n")
            fh.write("step,resname,resid,x_nm,y_nm,z_nm\n")

            groups, resid_map = _build_probe_groups(simulation.topology, allowed)
            anchor = _protein_anchor_indices(simulation.topology)

            # initial snapshot (step 0)
            st0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
            pos0 = st0.getPositions(asNumpy=False)
            ref_com = _com_of(pos0, anchor)
            _append_hotspots_csv(fh, 0, pos0, groups, resid_map, (0.0, 0.0, 0.0))
            fh.flush()

        if total_steps <= 0:
            return 0

        # segmented stepping
        remaining = int(total_steps)
        while remaining > 0:
            chunk = min(int(stride), remaining) if stride and stride > 0 else remaining
            simulation.step(chunk)
            step_counter += chunk
            remaining -= chunk

            if fh:
                st = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
                pos = st.getPositions(asNumpy=False)
                com_now = _com_of(pos, anchor) if anchor else (0.0, 0.0, 0.0)
                dx = (ref_com[0] - com_now[0]) if ref_com else 0.0
                dy = (ref_com[1] - com_now[1]) if ref_com else 0.0
                dz = (ref_com[2] - com_now[2]) if ref_com else 0.0
                _append_hotspots_csv(fh, step_counter, pos, groups, resid_map, (dx, dy, dz))
                fh.flush()

        return step_counter
    finally:
        if fh:
            try: fh.close()
            except Exception: pass

# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(
    description="Run OpenMM MD (RNA, protein, ligand(s), MixMD-from-Packmol)."
)
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--n-steps", type=int, default=1_000_000, help="Production MD steps (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed")

# Multi‑ligand convenience
parser.add_argument("--ligand-pattern", type=str, default="ligand_*.sdf",
                    help="Glob to find multiple ligands (default ligand_*.sdf). Falls back to ligand.sdf.")

# MixMD-from-Packmol options
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read placements CSV (build/<receptor>_mixmd_placements.csv) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACM,PHO,HAC,ACEA,PHOL,ACOH",
    help="Comma-separated probe residue names to import.")
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
    help="MD cubic box edge length; probes outside this cube are dropped.")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering Packmol centroids (nm).")
parser.add_argument("--hotspot-stride", type=int, default=1000,
    help="Stride (steps) for MixMD hotspot logging (segmented stepping).")
parser.add_argument("--hotspot-csv", type=str, default=None,
    help="Output CSV for MixMD hotspots (default: mixmd_hotspots_no_ligand.csv).")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo (ignore ligands)
if args.mixmd_from_packmol:
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---------------------------
# Load receptor
# ---------------------------
receptor_pdb = PDBFile(args.input_receptor)

# ---------------------------
# Force fields
# ---------------------------
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ---------------------------
# Build modeller & merge ligands or probes
# ---------------------------
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

lig_count = 0
if not args.no_ligand:
    # discover ligands
    lig_paths = sorted(Path(".").glob(args.ligand_pattern),
                       key=lambda p: (len(p.stem), p.stem))
    if not lig_paths and Path("ligand.sdf").exists():
        lig_paths = [Path("ligand.sdf")]
    if lig_paths:
        lig_count = _merge_ligands_into_modeller(modeller, forcefield, lig_paths)
        print(f"✅ Merged {lig_count} ligand(s): {[p.name for p in lig_paths]}")
    else:
        print("ℹ️  No ligand SDF found; running apo.")

# MixMD cosolvents (before solvation)
if args.mixmd_from_packmol:
    resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    added = _add_probes_from_packmol(
        modeller=modeller,
        forcefield=forcefield,
        receptor_path=receptor_path,
        sim_box_nm=float(args.mixmd_box_size_nm),
        edge_margin_nm=float(args.mixmd_edge_margin_nm),
        placements_csv=(Path(args.mixmd_placements_csv) if args.mixmd_placements_csv else None),
        resname_list=resnames
    )
    print(f"🧪 MixMD: added {added} probe instances.")

# Hydrogens after all non‑water molecules are present
modeller.addHydrogens(forcefield)

# Pre‑solvation snapshot
pre_path = _write_structure_auto(modeller.topology, modeller.positions,
                                 f"combined_receptor_ligand{suffix}")
print(f"✅ System ready for solvation → {pre_path}")

# ---------------------------
# Solvate in the requested MD box
# ---------------------------
box_edge = float(args.mixmd_box_size_nm if args.mixmd_from_packmol else 7.0)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(box_edge, box_edge, box_edge) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)
solv_path = _write_structure_auto(modeller.topology, modeller.positions,
                                  f"solvated_receptor_ligand{suffix}")
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm → {solv_path}")

# ---------------------------
# Create System
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

# Standard reporters (PDB is small/frequent; mmCIF is for whole snapshots we already wrote)
simulation.reporters.append(PDBReporter(f"minimized{suffix}.pdb", 100))
simulation.reporters.append(PDBReporter(f"nvt_equilibrated{suffix}.pdb", 100))
simulation.reporters.append(PDBReporter(f"npt_equilibrated{suffix}.pdb", 500))
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# ---------------------------
# Run (segmented with explicit hotspot logging)
# ---------------------------
print("🔹 Energy Minimization...")
simulation.minimizeEnergy()

print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * kelvin, seed)

# total segmentation stride
stride = max(1, int(args.hotspot_stride))

# MixMD hotspot CSV path (only if MixMD enabled)
hotspot_csv = None
probe_hint = None
if args.mixmd_from_packmol:
    probe_hint = {s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()}
    hotspot_csv = Path(args.hotspot_csv or f"mixmd_hotspots{suffix}.csv")

# NVT (500 steps)
_ = _run_segmented(simulation, total_steps=500, stride=stride,
                   csv_path=hotspot_csv, probe_names=probe_hint)

print("🔹 NPT Equilibration (5 ps)...")
_ = _run_segmented(simulation, total_steps=2500, stride=stride,
                   csv_path=hotspot_csv, probe_names=probe_hint)

print(f"🔥 Production MD: {args.n_steps:,} steps")
_ = _run_segmented(simulation, total_steps=int(args.n_steps), stride=stride,
                   csv_path=hotspot_csv, probe_names=probe_hint)

# Final snapshot
final_path = _write_structure_auto(simulation.topology,
                                   simulation.context.getState(getPositions=True).getPositions(),
                                   f"final_structure{suffix}")
print(f"🎉 MD complete → {final_path}")
