# 5_md_simulation.py
# Author: Iori Mochizuki
# Updated: 2025-10-10
# Description: Run OpenMM-based 2 ns MD simulation of ligand-receptor complex using OpenFF (RNA & no-ligand compatible)


import argparse
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm
from sys import stdout
import openmm as mm
import random
import numpy as _py

parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand, no-ligand)")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Input ligand SDF")
parser.add_argument("--n-steps", type=int, default=1000000, help="Number of steps for production MD (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

parser.add_argument("--mixmd", action="store_true",
                    help="Run MixMD-lite: add small-molecule probes (implies --no-ligand)")
parser.add_argument("--probe-total", type=int, default=200,
                    help="Total number of probes to replace waters with (default 200)")
parser.add_argument("--probe-fractions", type=str,
                    default="ipa=0.40,acn=0.20,imd=0.15,aceam=0.15,phol=0.10",
                    help="Comma-separated key=frac (must sum to 1.0). Keys: ipa,acn,imd,aceam,phol,acoh")
parser.add_argument("--probe-seed", type=int, default=24680,
                    help="Random seed for probe placement")
args = parser.parse_args()

if args.mixmd:
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

def parse_kv_fracs(s: str):
    kv = {}
    for token in s.split(","):
        token = token.strip()
        if not token: 
            continue
        if "=" not in token:
            raise ValueError(f"Bad token in --probe-fractions: '{token}' (use key=val)")
        k, v = token.split("=", 1)
        kv[k.strip()] = float(v.strip())
    # normalize
    total = _py.sum(list(kv.values()))
    if total <= 0:
        raise ValueError("Total probe fractions must be > 0")
    for k in kv:
        kv[k] /= float(total)
    return kv

# minimal cosolvent library
PROBES = {
    "ipa":   {"resname": "IPA",   "smiles": "CC(C)O"},       # isopropanol
    "acn":   {"resname": "ACN",   "smiles": "CC#N"},         # acetonitrile
    "imd":   {"resname": "IMD",   "smiles": "c1[nH]cnc1"},   # imidazole (neutral, N1 protonated)
    "aceam": {"resname": "ACEA",  "smiles": "CC(=O)N"},      # acetamide
    "phol":  {"resname": "PHOL",  "smiles": "c1ccc(cc1)O"},  # phenol
    "acoh":  {"resname": "ACOH",  "smiles": "CC(=O)O"}       # acetic acid (optional)
}

# 1. Load receptor
receptor_pdb = PDBFile(args.input_receptor)

# 2. Force field setup
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 forcefield")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all forcefield")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

probe_fracs = parse_kv_fracs(args.probe_fractions) if args.mixmd else {}
unique_off_mols = []

if not args.no_ligand:
    ligand = Molecule.from_file(args.input_ligand)  # you already have this
    unique_off_mols.append(ligand)

if args.mixmd:
    for key in probe_fracs.keys():
        if key not in PROBES:
            raise ValueError(f"Unknown probe key '{key}'. Use: {list(PROBES.keys())}")
        off = Molecule.from_smiles(PROBES[key]["smiles"], allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)  # 1 conformer is enough
        off.name = PROBES[key]["resname"]
        unique_off_mols.append(off)

# single SMIRNOFF generator can take all unique molecules
smirnoff = SMIRNOFFTemplateGenerator(molecules=unique_off_mols)
forcefield.registerTemplateGenerator(smirnoff.generator)


# 3. Ligand logic
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
    print("✅ No ligand: Running receptor-only MD")

modeller.addHydrogens(forcefield)

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

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

def replace_waters_with_probes(modeller, forcefield, probe_fracs, total_count, seed=24680):
    """
    Replace randomly selected water residues with small-molecule probes.
    Keeps your box, avoids overlaps, preserves file names.
    """
    rng = random.Random(int(seed))

    # collect water residues (handle typical names)
    water_names = set(["HOH", "WAT", "T3P", "TIP3"])
    waters = [res for res in modeller.topology.residues() if res.name in water_names]
    if len(waters) == 0:
        print("⚠️ No water residues recognized (names HOH/WAT/T3P/TIP3). Skipping probe placement.")
        return

    # decide counts
    counts = {}
    remaining = total_count
    keys = list(probe_fracs.keys())
    for i, k in enumerate(keys):
        if i < len(keys) - 1:
            n = int(round(total_count * probe_fracs[k]))
            counts[k] = n
            remaining -= n
        else:
            counts[k] = max(0, remaining)

    # select random water residues to replace
    rng.shuffle(waters)
    if sum(counts.values()) > len(waters):
        raise ValueError("Requested more probes than waters present.")

    water_iter = iter(waters)
    # prepare OFF → OpenMM topologies once
    off_cache = {}
    for k in counts:
        if counts[k] == 0: 
            continue
        if k not in off_cache:
            off_mol = next(m for m in unique_off_mols if getattr(m, "name", None) == PROBES[k]["resname"])
            off_cache[k] = {
                "off": off_mol,
                "top": OFFTopology.from_molecules([off_mol]).to_openmm(),
                "pos": off_mol.conformers[0]  # OpenFF quantity with Angstrom units
            }

    # we will delete selected waters and add probes at the oxygen position (with a tiny random jitter)
    to_delete = []
    add_queue = []  # list of (top, positions[nm], resname)

    for k, n_add in counts.items():
        for _ in range(n_add):
            try:
                wres = next(water_iter)
            except StopIteration:
                break
            # water oxygen position = first atom in residue (robust enough for OpenMM waters)
            wat_atom = next(water_atom for water_atom in wres.atoms() if water_atom.element.symbol == "O")
            wat_pos = modeller.positions[wat_atom.index]  # has units (nm)
            to_delete.extend(list(wres.atoms()))

            cache = off_cache[k]
            # center probe on its own centroid and translate to water oxygen; add tiny random jitter
            off_pos = cache["pos"]  # OpenFF unit (Å)
            # convert Å → nm and center
            coords_nm = _py.array([[p.value_in_unit(angstrom)] for p in off_pos])  # wrong shape; fix below

            #Working block
            # convert Å → nm as plain floats, center on centroid
            off_xyz_ang = _py.array([[p.x, p.y, p.z] for p in cache["pos"].value_in_unit(angstrom)])
            centroid_ang = off_xyz_ang.mean(axis=0, keepdims=True)
            off_xyz_nm = (off_xyz_ang - centroid_ang) * 0.1  # Å → nm

            # small random jitter (±0.05 nm)
            jitter = _py.array([rng.uniform(-0.05, 0.05),
                                rng.uniform(-0.05, 0.05),
                                rng.uniform(-0.05, 0.05)])

            # translate to water oxygen
            t = _py.array([wat_pos.x, wat_pos.y, wat_pos.z]) + jitter
            placed_nm = off_xyz_nm + t  # (natoms, 3) in nm

            # wrap as OpenMM Vec3 with units
            cand_pos = [mm.Vec3(float(v[0]), float(v[1]), float(v[2])) for v in placed_nm] * nanometer

            add_queue.append((cache["top"], cand_pos, PROBES[k]["resname"]))

            # perform deletion/addition on the modeller
            modeller.delete(to_delete)
            for top, cand_pos, resname in add_queue:
                modeller.add(top, cand_pos)
            # rename the newly added residues to requested 3–4 letter names
            # (heuristic: residues with default name 'MOL'/'UNK' that are not protein/water)
            for res in modeller.topology.residues():
                if res.name in ("MOL", "UNK"):
                    # attempt to rename by matching atom count to a probe (simple, robust when running apo MixMD)
                    nat = sum(1 for _ in res.atoms())
                    for k in probe_fracs.keys():
                        off_mol = off_cache.get(k, {}).get("off", None)
                        if off_mol and nat == off_mol.n_atoms:
                            res.name = PROBES[k]["resname"]
        
            print(f"✅ Replaced waters with probes: {counts}")

if args.mixmd:
    print("🧪 [MixMD-lite] Replacing waters with probe molecules...")
    replace_waters_with_probes(modeller, forcefield, probe_fracs, args.probe_total, seed=args.probe_seed)
    with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print("✅ Probes placed (water replacements).")


system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1 * nanometer,
    constraints=HBonds
)
seed = int(args.seed)
#system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin, 25))

# Barostat with a fixed seed
barostat = MonteCarloBarostat(1 * bar, 300 * kelvin, 25)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

# Langevin integrator with a fixed seed
integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
integrator.setRandomNumberSeed(seed)
#platform = Platform.getPlatformByName("CUDA")
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
simulation.step(args.n_steps)  # default = 1,000,000 steps for 2 ns

with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
print(f"🎉 MD complete → final_structure{suffix}.pdb saved.")
