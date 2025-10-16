# 5_md_simulation.py
# Author: Iori Mochizuki
# Updated: 2025-10-15
# Description: OpenMM 2 ns MD (protein/RNA/ligand) with optional MixMD-from-Packmol import.
# - Keeps original behavior and file names.
# - When --mixmd-from-packmol is supplied:
#     * implies --no-ligand
#     * reads build/<receptor_stem>_mixmd.pdb and ..._mixmd_placements.csv (or user-provided paths)
#     * filters probe placements to a user-sized MD box (default 7.0 nm)
#     * parametrizes probes via OpenFF/SMIRNOFF and adds them before solvation.

import argparse, csv
from pathlib import Path
from sys import stdout

#import numpy as _np
import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# --- MixMD-from-Packmol helpers (pure OpenMM, no NumPy) -----------------------

PROBE_LIBRARY = {
    "IPA":  "CC(C)O",          # isopropanol
    "ACN":  "CC#N",            # acetonitrile
    "IMD":  "c1[nH]cnc1",      # imidazole (neutral N1-H)
    "ACEA": "CC(=O)N",         # acetamide
    "PHOL": "c1ccc(cc1)O",     # phenol
    "ACOH": "CC(=O)O",         # acetic acid
}

def _default_packmol_paths(receptor_path: Path):
    root = f"{receptor_path.stem}_mixmd"
    packmol_pdb = Path("build") / f"{root}.pdb"
    placements  = Path("build") / f"{root}_placements.csv"
    return packmol_pdb, placements

def _openff_conf_to_angstrom_list(off_mol):
    """Return a list of (x,y,z) floats in Å from an OpenFF Molecule conformer using plain Python loops."""
    conf = off_mol.conformers[0]
    # Try OpenMM-unit path first
    try:
        arr = conf.value_in_unit(angstrom)  # N×3
        out = []
        for row in arr:
            out.append((float(row[0]), float(row[1]), float(row[2])))
        return out
    except Exception:
        # Fallback: assume nested sequences
        out = []
        for row in conf:
            out.append((float(row[0]), float(row[1]), float(row[2])))
        return out

def _positions_for_centroid_nm(off_mol, centroid_nm_tuple):
    """
    Build OpenMM positions (Quantity of list[Vec3] in nanometer) for `off_mol`
    centered on centroid_nm_tuple (x,y,z in nm). Pure Python, no NumPy.
    """
    coords_A = _openff_conf_to_angstrom_list(off_mol)   # [(x,y,z)_Å]
    n = len(coords_A)
    cx = sum(p[0] for p in coords_A) / n
    cy = sum(p[1] for p in coords_A) / n
    cz = sum(p[2] for p in coords_A) / n
    tx, ty, tz = (float(centroid_nm_tuple[0]),
                  float(centroid_nm_tuple[1]),
                  float(centroid_nm_tuple[2]))
    placed = []
    for (xA, yA, zA) in coords_A:
        xnm = (xA - cx) * 0.1 + tx   # Å→nm and translate
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
    - `placements_csv`: optional override; default is build/<receptor>_mixmd_placements.csv
    - `resname_list`: optional whitelist of probe residue names (e.g., ["IPA","ACN"])
    Returns: number of probes actually added.
    """
    # Resolve default paths (we only need the CSV to place centroids)
    if placements_csv is None:
        _, default_csv = _default_packmol_paths(receptor_path)
        placements_csv = default_csv

    print(f"ℹ️  MixMD inputs: CSV={placements_csv}")

    if not placements_csv.exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")

    # Half-box for filtering
    half = 0.5 * float(sim_box_nm)
    keep_lo = - (half - float(edge_margin_nm))
    keep_hi = + (half - float(edge_margin_nm))

    # Optional resname whitelist
    whitelist = None
    if resname_list:
        whitelist = {r.upper().strip() for r in resname_list}

    # Parse rows and crop to the MD box
    rows = []
    for rec in _read_centroids_csv(placements_csv):
        resname = rec["resname"]
        if resname not in PROBE_LIBRARY:
            continue
        if whitelist is not None and resname not in whitelist:
            continue
        x = rec["x_nm"]; y = rec["y_nm"]; z = rec["z_nm"]
        if (keep_lo <= x <= keep_hi) and (keep_lo <= y <= keep_hi) and (keep_lo <= z <= keep_hi):
            rows.append((resname, (x, y, z)))

    if not rows:
        print("⚠️  No probe centroids inside the target box; continuing with receptor-only.")
        return 0

    # Prepare OFF molecules and register SMIRNOFF once
    used_res = sorted({res for (res, _) in rows})
    off_mols = []
    for res in used_res:
        smi = PROBE_LIBRARY[res]
        off = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = res  # readability only
        off_mols.append(off)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    # Place probes
    added = 0
    for (resname, cen_nm) in rows:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)  # list of Vec3*nm
        modeller.add(top, pos)
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
                    help="MD cubic box edge length; probes outside this cube are dropped (default 7.0 nm).")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACEA,PHOL,ACOH",
                    help="Comma-separated probe residue names to import (default: all).")
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read Packmol mixture (build/<receptor>_mixmd*.{pdb,csv}) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-packmol-pdb", type=str, default=None,
    help="Override path to Packmol mixture PDB (default: build/<receptor>_mixmd.pdb).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV written by 3c (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering centroids (default 0.15 nm).")


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

# We may add SMIRNOFF later if probes (or a ligand) are present.
unique_off_mols = []

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

# If Packmol probes are requested, add them now (before solvation)
if args.mixmd_from_packmol:
    resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    added = _add_probes_from_packmol(
        modeller=modeller,
        forcefield=forcefield,
        receptor_path=receptor_path,
        packmol_pdb=Path(args.mixmd_packmol_pdb) if args.mixmd_packmol_pdb else None,
        placements_csv=Path(args.mixmd_placements_csv) if args.mixmd_placements_csv else None,
        box_size_nm=float(args.mixmd_box_size_nm),
        resname_list=resnames
    )
    print(f"🧪 MixMD: added {added} probe instances.")

# Hydrogens after all non‑water molecules are present
modeller.addHydrogens(forcefield)

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# If requested, add Packmol probes now (before solvation)
if getattr(args, "mixmd_from_packmol", False):
    print("🧪 [MixMD-from-Packmol] Using placements CSV to add probes (apo).")
    receptor_path = Path(args.input_receptor)
    # Use the same box length you will pass to addSolvent (7.0 nm in your script)
    SIM_BOX_NM = 7.0
    _add_probes_from_packmol(
        modeller,
        forcefield,
        receptor_path,
        sim_box_nm=SIM_BOX_NM,
        edge_margin_nm=float(getattr(args, "mixmd_edge_margin_nm", 0.15))
    )


# ---------------------------
# Solvate in the requested MD box (default 7.0 nm)
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
