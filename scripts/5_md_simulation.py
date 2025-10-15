# 5_md_simulation.py
# Author: Iori Mochizuki (pipeline owner) + collaborator patch
# Updated: 2025-10-15
# Description: OpenMM MD (protein/RNA ± ligand). Optional MixMD-from-Packmol mode.
# - Defaults unchanged from your original script.
# - --mixmd-from-packmol reads build/<receptor_stem>_mixmd.{pdb,csv} from 3c and places probes.

import argparse, os, sys, csv
from pathlib import Path

from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
import openmm as mm

from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---- CLI --------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand, no-ligand)")

# original flags (unchanged)
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Input ligand SDF")
parser.add_argument("--n-steps", type=int, default=1000000, help="Production steps (default: 1,000,000 ≈ 2 ns)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# new: Packmol/MixMD mode
parser.add_argument("--mixmd-from-packmol", action="store_true",
                    help="Add Packmol probes before solvation (implies --no-ligand by default).")
parser.add_argument("--packmol-prefix", type=str, default=None,
                    help="Prefix to Packmol outputs (default: build/<receptor_stem>_mixmd)")
parser.add_argument("--probe-smiles", type=str, default="ipa=CC(C)O,acn=CC#N,imd=c1[nH]cnc1,aceam=CC(=O)N,phol=c1ccc(cc1)O,acoh=CC(=O)O",
                    help="Comma list key=SMILES for probes (keys must match resnames in 3c: IPA,ACN,IMD,ACEA,PHOL,ACOH).")

args = parser.parse_args()

# In MixMD mode, we run apo unless the user explicitly wants a ligand
if args.mixmd_from_packmol:
    if not args.no_ligand:
        print("ℹ️  --mixmd-from-packmol requested → running APO (imply --no-ligand).")
        args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---- Helpers ----------------------------------------------------------------
def _default_packmol_paths(receptor_path: Path, tag: str = "mixmd"):
    """
    Given receptor_cleaned.pdb -> defaults to:
      build/receptor_cleaned_mixmd.pdb
      build/receptor_cleaned_mixmd_placements.csv
    """
    build = Path("build")
    prefix = build / f"{Path(receptor_path).stem}_{tag}"

    packmol_pdb = prefix.with_suffix(".pdb")                     # OK: real suffix
    placements_csv = prefix.with_name(prefix.name + "_placements.csv")  # tail via with_name
    return packmol_pdb, placements_csv

# Resolve defaults if flags are not set
if args.mixmd_packmol_pdb is None or args.mixmd_placements_csv is None:
    packmol_pdb, packmol_csv = _default_packmol_paths(receptor_path)
else:
    packmol_pdb = Path(args.mixmd_packmol_pdb)
    packmol_csv = Path(args.mixmd_placements_csv)

# Helpful diagnostics
print(f"ℹ️  MixMD inputs: PDB={packmol_pdb} | CSV={packmol_csv}")

if not packmol_pdb.exists() or not packmol_csv.exists():
    raise FileNotFoundError(
        "Packmol outputs not found.\n"
        f"  Expected PDB: {packmol_pdb}\n"
        f"  Expected CSV: {packmol_csv}\n"
        "Tip: 3c writes into the ./build/ folder by default. Run 3c in the SAME project "
        "directory as 5, or pass explicit paths with:\n"
        "  --mixmd-packmol-pdb build/<stem>_mixmd.pdb "
        "--mixmd-placements-csv build/<stem>_mixmd_placements.csv"
    )


def _parse_cryst1_box_nm(pdb_path: Path):
    """Return cubic box length (nm) from CRYST1 if present; else None."""
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if line.startswith("CRYST1") and len(line) >= 54:
                    # CRYST1 a b c alpha beta gamma
                    aA = float(line[6:15]); bA = float(line[15:24]); cA = float(line[24:33])
                    # Use the max edge to ensure everything fits; convert Å→nm
                    L_nm = max(aA, bA, cA) * 0.1
                    return L_nm
    except Exception:
        pass
    return None

def _load_packmol_centroids(csv_path: Path):
    """Read centroids CSV produced by 3c (resname,resid,x_nm,y_nm,z_nm)."""
    centroids = []  # list of dicts
    with open(csv_path, "r", newline="") as fp:
        rdr = csv.DictReader(fp)
        for row in rdr:
            centroids.append({
                "resname": row["resname"].strip().upper(),
                "resid": int(row["resid"]),
                "x": float(row["x_nm"]),
                "y": float(row["y_nm"]),
                "z": float(row["z_nm"]),
            })
    if not centroids:
        raise RuntimeError(f"No centroids parsed from {csv_path}")
    return centroids

def _parse_probe_smiles_map(s: str):
    """
    'ipa=CC(C)O,acn=CC#N,...' → dict:
      {'IPA': 'CC(C)O', 'ACN': 'CC#N', ... }
    Resnames we use in 3c: IPA, ACN, IMD, ACEA, PHOL, ACOH
    """
    m = {}
    for tok in s.split(","):
        if not tok.strip():
            continue
        k, smi = tok.split("=", 1)
        k = k.strip().lower()
        smi = smi.strip()
        # map key -> resname used by 3c
        key_to_res = {
            "ipa": "IPA", "acn": "ACN", "imd": "IMD",
            "aceam": "ACEA", "phol": "PHOL", "acoh": "ACOH"
        }
        if k not in key_to_res:
            raise ValueError(f"Unknown probe key '{k}' in --probe-smiles")
        m[key_to_res[k]] = smi
    return m

def _build_probe_prototypes(resname_to_smiles: dict):
    """
    For each resname used in Packmol, create an OpenFF Molecule with 1 conformer,
    centered at origin. Return dict resname -> dict(off, top, centered_nm).
    """
    out = {}
    for resname, smi in resname_to_smiles.items():
        off = Molecule.from_smiles(smi, allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = resname
        top = OFFTopology.from_molecules([off]).to_openmm()
        # center conformer to its centroid (Å→nm)
        conf = off.conformers[0]  # Quantities in Å
        coordsA = conf.value_in_unit(angstrom)  # ndarray (n,3)
        cx, cy, cz = coordsA.mean(axis=0)
        centered_nm = (coordsA - [cx, cy, cz]) * 0.1  # nm
        out[resname] = {"off": off, "top": top, "centered_nm": centered_nm}
    return out

def _add_packmol_probes_to_modeller(modeller: Modeller,
                                    prototypes: dict,
                                    centroids: list):
    """
    Add one instance per centroid: translate prototype coords to centroid (nm),
    convert to Vec3 (nm), then modeller.add().
    """
    # Build once: SMIRNOFFTemplateGenerator will parametrize any residues
    # matching these molecules at createSystem().
    # (Registration is done in the main body, before createSystem.)
    # Here we only add positions/topologies.
    added = 0
    for c in centroids:
        resname = c["resname"]
        if resname not in prototypes:
            # skip centroids whose resname was not requested in --probe-smiles
            continue
        proto = prototypes[resname]
        base = proto["centered_nm"]  # (n,3) nm
        tx, ty, tz = float(c["x"]), float(c["y"]), float(c["z"])
        placed = [mm.Vec3(float(x + tx), float(y + ty), float(z + tz)) for x, y, z in base]
        modeller.add(proto["top"], placed * nanometer)
        added += 1
    print(f"✅ Added {added} Packmol probe instances to modeller.")

# ---- Load receptor & forcefield ---------------------------------------------
receptor_path = Path(args.input_receptor)
receptor_pdb = PDBFile(receptor_path.as_posix())

if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ---- Optional ligand (unchanged behavior) -----------------------------------
unique_off_mols = []
if not args.no_ligand:
    ligand = Molecule.from_file(args.input_ligand)
    unique_off_mols.append(ligand)

# ---- Optional Packmol probes (new) ------------------------------------------
packmol_pdb = None
packmol_csv = None
box_from_packmol_nm = None
prototypes = None

if args.mixmd_from_packmol:
    # Resolve default Packmol outputs for this receptor
    if args.packmol_prefix is None:
        packmol_pdb, packmol_csv = _default_packmol_paths(receptor_path)
    else:
        pfx = Path(args.packmol_prefix)
        packmol_pdb = pfx.with_suffix(".pdb")
        # placements CSV name follows our 3c script convention
        if pfx.suffix:
            # if user passed already ending with .pdb or .inp, derive stem
            stem = pfx.with_suffix("").name
            packmol_csv = pfx.with_suffix("").with_name(stem + "_placements.csv")
        else:
            packmol_csv = pfx.with_name(pfx.name + "_placements.csv")

    if not packmol_pdb.exists():
        raise FileNotFoundError(f"Packmol PDB not found: {packmol_pdb}")
    if not packmol_csv.exists():
        raise FileNotFoundError(f"Packmol placements CSV not found: {packmol_csv}")

    # parse box
    box_from_packmol_nm = _parse_cryst1_box_nm(packmol_pdb)
    if box_from_packmol_nm is None:
        print("⚠️  No CRYST1 found in Packmol PDB; will fall back to 7.0 nm box.")

    # which probes to parameterize
    res_to_smi = _parse_probe_smiles_map(args.probe_smiles)
    prototypes = _build_probe_prototypes(res_to_smi)
    unique_off_mols.extend([p["off"] for p in prototypes.values()])

# Register a single SMIRNOFF generator that knows about all nonstandard molecules
if unique_off_mols:
    smirnoff = SMIRNOFFTemplateGenerator(molecules=unique_off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

# ---- Build modeller ----------------------------------------------------------
if not args.no_ligand:
    # Holo: add receptor + ligand as in your original script
    ligand_positions = to_openmm(ligand.conformers[0])
    ligand_top = OFFTopology.from_molecules([ligand]).to_openmm()
    modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)
    modeller.add(ligand_top, ligand_positions)
    print("✅ Ligand merged")
else:
    modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)
    print("✅ Apo: receptor only (so far)")

# NEW: add Packmol probes before solvation
if args.mixmd_from_packmol:
    centroids = _load_packmol_centroids(packmol_csv)
    _add_packmol_probes_to_modeller(modeller, prototypes, centroids)

# Hydrogens (proteins/RNA; probes already have explicit H from OpenFF)
modeller.addHydrogens(forcefield)

# Write pre-solvation combined structure (same filename convention)
with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# Solvate (box from Packmol if available; else 7.0 nm as before)
if box_from_packmol_nm is not None:
    L = float(box_from_packmol_nm)
    box_vec = (L, L, L) * nanometer
else:
    box_vec = (7.0, 7.0, 7.0) * nanometer

modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=box_vec,
    ionicStrength=0.15 * molar,
    neutralize=True
)
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("✅ Solvated system ready.")

# ---- System & Simulation (unchanged) ----------------------------------------
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1 * nanometer,
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
simulation.reporters.append(StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, temperature=True))
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
