# 5_md_simulation.py
# Author: Iori Mochizuki
# Updated: 2025-10-16
# Description: OpenMM 2 ns MD (protein/RNA/ligand) with optional MixMD-from-Packmol import.
# - When --mixmd-from-packmol is supplied:
#     * implies --no-ligand
#     * reads build/<receptor_stem>_mixmd_placements.csv (or user-provided path)
#     * filters probe placements to a user-sized MD box (default 7.0–8.0 nm)
#     * parametrizes probes via OpenFF/SMIRNOFF and adds them before solvation
#     * RENAMES added residues to probe resnames (IPA/ACN/IMD/ACEA/PHOL/ACOH)
#     * logs probe centroids during MD to a CSV (hotspot reporter)

import argparse, csv
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
# Hotspot reporter
# ---------------------------
class ProbeHotspotReporter:
    """
    Write probe residue centroids (nm) every 'reportInterval' steps to CSV.
    Columns: step,resname,resid,x_nm,y_nm,z_nm

    Auto-detects probe residues present in the provided topology by intersecting
    residue names with PROBE_NAME_SET (accepts 3-letter canonical and 4-letter aliases).
    """
    def __init__(self, file_path, reportInterval, topology, probe_resnames=None):
        self.file = open(file_path, "w")
        self.reportInterval = int(reportInterval)

        # Build the set of names to track: auto-discover ∩ optional user list
        present = { (res.name or "").upper().strip() for res in topology.residues() }
        auto_names = present & PROBE_NAME_SET

        if probe_resnames:
            user_names = { s.upper().strip() for s in probe_resnames }
            self.probe_res = auto_names | (user_names & PROBE_NAME_SET)
        else:
            self.probe_res = auto_names

        # Build index groups once
        groups = []
        resid_map = []  # (resname, resid)
        for res in topology.residues():
            rn = (res.name or "").upper().strip()
            if rn in self.probe_res:
                idxs = [a.index for a in res.atoms()]
                if idxs:
                    groups.append(idxs)
                    # PDB-like residue id if available; otherwise 0
                    resid_map.append((rn, int(res.id) if res.id is not None else 0))

        self.groups = groups
        self.resid_map = resid_map

        # header + quick diagnostic in the CSV (as a comment)
        self.file.write("# tracking_resnames=" + ";".join(sorted(self.probe_res)) + "\n")
        self.file.write("step,resname,resid,x_nm,y_nm,z_nm\n")
        self.file.flush()

    def describeNextReport(self, simulation):
        # (steps, positions, velocities, forces, energies, enforcePeriodicBox)
        return (self.reportInterval, True, False, False, False, True)

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=False)  # list of Vec3 (nm)

        # --- robust step retrieval across OpenMM versions ---
        try:
            step = int(simulation.currentStep)      # OpenMM ≥ 7.5+
        except Exception:
            # Fallback: derive from time/stepsize with units
            try:
                t  = state.getTime()                           # has units
                dt = simulation.integrator.getStepSize()       # has units
                step = int(round((t / dt)))
            except Exception:
                step = -1  # last resort; still log coordinates
        # -----------------------------------------------------

        for (idxs, (resname, resid)) in zip(self.groups, self.resid_map):
            sx = sy = sz = 0.0
            n = 0
            for i in idxs:
                v = pos[i]  # Vec3 (nm units)
                sx += float(v.x); sy += float(v.y); sz += float(v.z)
                n += 1
            if n > 0:
                self.file.write(f"{step},{resname},{resid},{sx/n:.5f},{sy/n:.5f},{sz/n:.5f}\n")
        self.file.flush()

    def __del__(self):
        try: self.file.close()
        except Exception: pass


# --- MixMD-from-Packmol helpers (pure OpenMM, no NumPy) -----------------------
# Canonical 3-letter names used by 3c/Packmol, plus 4-letter aliases so 5 accepts both.
PROBE_LIBRARY = {
    # canonical 3-letter → SMILES
    "IPA":  "CC(C)O",           # isopropanol
    "ACN":  "CC#N",             # acetonitrile
    "IMD":  "c1[nH]cnc1",       # imidazole (neutral, N1–H)
    "ACM":  "CC(=O)N",          # acetamide
    "PHO":  "c1ccc(cc1)O",      # phenol
    "HAC":  "CC(=O)O",          # acetic acid
    # aliases (legacy 4-letter spellings that may appear in older runs/docs)
    "ACEA": "CC(=O)N",          # acetamide (alias of ACM)
    "PHOL": "c1ccc(cc1)O",      # phenol   (alias of PHO)
    "ACOH": "CC(=O)O",          # acetic acid (alias of HAC)
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())


def _default_packmol_paths(receptor_path: Path):
    root = f"{receptor_path.stem}_mixmd"
    placements  = Path("build") / f"{root}_placements.csv"
    return placements

def _openff_conf_to_angstrom_list(off_mol):
    """
    Return a list of (x,y,z) floats in Å from an OpenFF Molecule conformer.
    Pure Python loops; works with Pint-backed units.
    """
    conf = off_mol.conformers[0]  # OpenFF quantity w/ units
    # Preferred path: Pint → Å magnitudes
    try:
        arr = conf.to("angstrom").m  # N×3 numeric
        out = [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
        return out
    except Exception:
        # Fallback via OpenMM units
        q = to_openmm(conf)  # Quantity(list(Vec3), angstrom)
        out = []
        for v in q.value_in_unit(angstrom):
            out.append((float(v[0]), float(v[1]), float(v[2])))
        return out

def _positions_for_centroid_nm(off_mol, centroid_nm_tuple):
    """
    Build OpenMM positions (list[Vec3]*nanometer) for `off_mol`
    centered on centroid_nm_tuple (x,y,z in nm). Pure Python math.
    """
    coords_A = _openff_conf_to_angstrom_list(off_mol)   # list of (x,y,z) in Å
    n = len(coords_A)
    if n == 0:
        raise ValueError("OFF molecule has no coordinates.")

    # centroid in Å (plain floats)
    sx = sy = sz = 0.0
    for xA, yA, zA in coords_A:
        sx += xA; sy += yA; sz += zA
    cx = sx / n; cy = sy / n; cz = sz / n

    tx, ty, tz = (float(centroid_nm_tuple[0]),
                  float(centroid_nm_tuple[1]),
                  float(centroid_nm_tuple[2]))

    placed = []
    for xA, yA, zA in coords_A:
        xnm = (xA - cx) * 0.1 + tx   # Å→nm + translate
        ynm = (yA - cy) * 0.1 + ty
        znm = (zA - cz) * 0.1 + tz
        placed.append(mm.Vec3(xnm, ynm, znm))
    return [p * nanometer for p in placed]

def _read_centroids_csv(csv_path: Path):
    """Read placements CSV written by 3c; returns list of dicts."""
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

def _add_probes_from_packmol(modeller,
                             forcefield,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None = None,
                             resname_list: list[str] | None = None) -> int:
    """
    Add probes at Packmol centroids cropped to the intended MD box.
    Critically: rename the newly added residues to the probe resname
    so downstream reporters can find them.
    Returns: number of probes actually added.
    """
    # Resolve CSV path
    if placements_csv is None:
        placements_csv = _default_packmol_paths(receptor_path)

    print(f"ℹ️  MixMD inputs: CSV={placements_csv}")
    if not placements_csv.exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")

    # Half-box for filtering
    half = 0.5 * float(sim_box_nm)
    keep_lo = - (half - float(edge_margin_nm))
    keep_hi = + (half - float(edge_margin_nm))

    # Resname whitelist
    whitelist = {r.upper().strip() for r in resname_list} if resname_list else None

    # Collect rows inside the MD cube
    rows = []
    for rec in _read_centroids_csv(placements_csv):
        resname = rec["resname"]
        if resname not in PROBE_LIBRARY:
            continue
        if whitelist and resname not in whitelist:
            continue
        x, y, z = rec["x_nm"], rec["y_nm"], rec["z_nm"]
        if (keep_lo <= x <= keep_hi) and (keep_lo <= y <= keep_hi) and (keep_lo <= z <= keep_hi):
            rows.append((resname, rec["resid"], (x, y, z)))

    if not rows:
        print("⚠️  No probe centroids inside the target box; continuing with receptor-only.")
        return 0

    # Prepare OFF molecules and register SMIRNOFF once
    used_res = sorted({res for (res, _, _) in rows})
    off_mols = []
    for res in used_res:
        smi = PROBE_LIBRARY[res]
        off = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = res
        off_mols.append(off)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    # Place probes and RENAME newly added residues
    def _res_list():
        return list(modeller.topology.residues())

    added = 0
    for (resname, pack_resid, cen_nm) in rows:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)

        # remember residues before add
        before_res = len(_res_list())
        modeller.add(top, pos)
        after = _res_list()

        # newly added residues are the tail slice
        new_res = after[before_res:]
        for res in new_res:
            try:
                res.name = resname
            except Exception:
                pass
            # keep Packmol resid if possible (purely cosmetic / useful for logs)
            try:
                res.id = str(pack_resid)
            except Exception:
                pass

        added += 1

    print(f"✅ Placed {added} probe molecules (cropped to {sim_box_nm:.1f} nm box).")
    return added

# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand, no-ligand, MixMD-from-Packmol).")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Input ligand SDF")
parser.add_argument("--n-steps", type=int, default=1000000, help="Number of steps for production MD (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# MixMD-from-Packmol options
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
                    help="MD cubic box edge length; probes outside this cube are dropped.")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACEA,PHOL,ACOH",
                    help="Comma-separated probe residue names to import.")
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read placements CSV (build/<receptor>_mixmd_placements.csv) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering centroids.")
parser.add_argument("--hotspot-csv", type=str, default=None,
    help="If set with --mixmd-from-packmol, write probe centroids every --hotspot-stride to this CSV (default: mixmd_hotspots_no_ligand.csv).")
parser.add_argument("--hotspot-stride", type=int, default=1000,
    help="Stride (steps) for hotspot centroid logging.")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo
if getattr(args, "mixmd_from_packmol", False):
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
# Build modeller (ligand or apo)
# ---------------------------
if not args.no_ligand:
    ligand = Molecule.from_file(args.input_ligand)
    ligand_positions = to_openmm(ligand.conformers[0])
    ligand_top = OFFTopology.from_molecules([ligand]).to_openmm()
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[ligand])
    forcefield.registerTemplateGenerator(smirnoff.generator)
    modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)
    modeller.add(ligand_top, ligand_positions)
    print("✅ Ligand merged")
else:
    modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)
    print("✅ No ligand: Running receptor-only (apo) setup")

# ---------------------------
# If Packmol probes are requested, add them now (before solvation)
# ---------------------------
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

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

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
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm")

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

# --- Hotspot reporter (only if MixMD) ---
if args.mixmd_from_packmol:
    # Let the reporter auto-discover present probe residue names; still pass user list as hints.
    probe_resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    hotspot_csv = args.hotspot_csv or f"mixmd_hotspots{suffix}.csv"

    # Create reporter
    hrep = ProbeHotspotReporter(
        hotspot_csv,
        reportInterval=int(args.hotspot_stride),
        topology=simulation.topology,
        probe_resnames=probe_resnames
    )

    # Force an initial snapshot at step 0 so the CSV is never empty
    state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    hrep.report(simulation, state0)

    # Now add it to the reporting chain for the run
    simulation.reporters.append(hrep)


# ---------------------------
# Run
# ---------------------------
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

with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
print(f"🎉 MD complete → final_structure{suffix}.pdb saved.")
