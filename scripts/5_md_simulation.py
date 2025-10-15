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

import numpy as _np
import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---------------------------
# Probe library (resname -> SMILES)
# ---------------------------
PROBE_MAP = {
    "IPA":  "CC(C)O",          # isopropanol
    "ACN":  "CC#N",            # acetonitrile
    "IMD":  "c1[nH]cnc1",      # imidazole (neutral)
    "ACEA": "CC(=O)N",         # acetamide
    "PHOL": "c1ccc(cc1)O",     # phenol
    "ACOH": "CC(=O)O",         # acetic acid (optional)
}

# ---------------------------
# Helpers
# ---------------------------
def _default_packmol_paths(receptor_path: Path):
    """Return default Packmol output paths for a given receptor path."""
    stem = receptor_path.stem
    build = Path("build")
    return build / f"{stem}_mixmd.pdb", build / f"{stem}_mixmd_placements.csv"

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

def _make_off_molecules(resnames):
    """Build OFF Molecule objects (with 1 conformer) for the requested probe resnames."""
    off = {}
    for res in resnames:
        if res not in PROBE_MAP:
            raise ValueError(f"Unknown probe resname '{res}'. Allowed: {list(PROBE_MAP.keys())}")
        m = Molecule.from_smiles(PROBE_MAP[res], allow_undefined_stereo=True)
        m.generate_conformers(n_conformers=1)
        # keep a readable name; residue name in topology may still appear as MOL/UNL, which is fine
        m.name = res
        off[res] = m
    return off

def _positions_for_centroid(off_mol: Molecule, centroid_nm):
    """
    Create OpenMM positions for 'off_mol' centered at 'centroid_nm' (x,y,z in nm).
    Returns a unit.Quantity(list(Vec3), nanometer).
    """
    # OFF conformer coordinates come with units (Å) via to_openmm
    conf_angs = to_openmm(off_mol.conformers[0])  # Quantity(list(Vec3), angstrom)
    # Convert to plain float array in Å
    arr_A = _np.asarray([[v.x, v.y, v.z] for v in conf_angs.value_in_unit(angstrom)], dtype=float)
    # Center molecule at its own centroid
    arr_A -= arr_A.mean(axis=0, keepdims=True)
    # Convert Å -> nm
    arr_nm = arr_A * 0.1
    # Translate to centroid
    tx, ty, tz = centroid_nm
    arr_nm[:, 0] += tx
    arr_nm[:, 1] += ty
    arr_nm[:, 2] += tz
    # Wrap back to Vec3 and attach units
    vecs = [mm.Vec3(float(x), float(y), float(z)) for x, y, z in arr_nm]
    return Quantity(vecs, nanometer)

def _add_probes_from_packmol(modeller: Modeller,
                             forcefield: ForceField,
                             receptor_path: Path,
                             packmol_pdb: Path|None,
                             placements_csv: Path|None,
                             box_size_nm: float,
                             resname_list: list[str]):
    """
    Import Packmol probes but only keep those inside the target MD box (centered at 0, cube of edge box_size_nm).
    Adds SMIRNOFF parameters for all used probes to 'forcefield' and adds molecules to 'modeller'.
    """
    # Resolve defaults
    if packmol_pdb is None or placements_csv is None:
        default_pdb, default_csv = _default_packmol_paths(receptor_path)
        packmol_pdb = packmol_pdb or default_pdb
        placements_csv = placements_csv or default_csv

    if not Path(placements_csv).exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")
    if not Path(packmol_pdb).exists():
        raise FileNotFoundError(f"Packmol mixture PDB not found: {packmol_pdb}")

    # Read centroids
    rows = _read_centroids_csv(placements_csv)

    # Filter by resname and by MD box
    half = 0.5 * float(box_size_nm)
    wanted = set(r.upper() for r in resname_list)
    kept = [r for r in rows
            if (r["resname"] in wanted
                and (-half <= r["x_nm"] <= half)
                and (-half <= r["y_nm"] <= half)
                and (-half <= r["z_nm"] <= half))]

    if len(kept) == 0:
        print("⚠️  No probes fall inside the requested MD box; continuing without probes.")
        return 0

    # OFF molecules and SMIRNOFF registration
    unique_res = sorted(set(r["resname"] for r in kept))
    off_mols = _make_off_molecules(unique_res)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=list(off_mols.values()))
    forcefield.registerTemplateGenerator(smirnoff.generator)

    # Add each instance
    added = 0
    for r in kept:
        res = r["resname"]
        off = off_mols[res]
        top = OFFTopology.from_molecules([off]).to_openmm()
        pos = _positions_for_centroid(off, (r["x_nm"], r["y_nm"], r["z_nm"]))
        modeller.add(top, pos)
        added += 1

    print(f"✅ Imported {added} probe instances from Packmol (kept inside ±{half:.2f} nm).")
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
parser.add_argument("--mixmd-from-packmol", action="store_true",
                    help="Import cosolvent probes from Packmol (3c output). Implies --no-ligand.")
parser.add_argument("--mixmd-packmol-pdb", type=str, default=None,
                    help="Path to Packmol mixture PDB (default: build/<stem>_mixmd.pdb).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
                    help="Path to placements CSV (default: build/<stem>_mixmd_placements.csv).")
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
                    help="MD cubic box edge length; probes outside this cube are dropped (default 7.0 nm).")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACEA,PHOL,ACOH",
                    help="Comma-separated probe residue names to import (default: all).")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo
if args.mixmd_from_packmol:
    print("ℹ️  --mixmd-from-packmol requested → running APO (implies --no-ligand).")
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
