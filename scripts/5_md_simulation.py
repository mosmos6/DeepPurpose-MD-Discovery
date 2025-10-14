# 5_md_simulation.py (add-only, non-breaking update)
# Author: Iori Mochizuki (base) / minor hook by ChatGPT
# Updated: 2025-10-14
# Description: Original 2 ns OpenMM MD with an optional --packmol-probes-pdb hook

import argparse
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm
from sys import stdout
import openmm as mm

# ---- NEW: optional Packmol hook flags (universal) ----------------------------
PROBE_DEF = {
    # resname : (smiles)
    "IPA":  "CC(C)O",
    "ACN":  "CC#N",
    "IMD":  "c1[nH]cnc1",
    "ACEA": "CC(=O)N",
    "PHOL": "c1ccc(cc1)O",
    "ACOH": "CC(=O)O",
}
# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, ligand, no-ligand)")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Input ligand SDF")
parser.add_argument("--n-steps", type=int, default=1000000, help="Number of steps for production MD (default: 1,000,000)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# ---- NEW: if present, add probes (no changes if omitted) ---------------------
parser.add_argument("--packmol-probes-pdb", type=str, default=None,
                    help="Optional Packmol mixture PDB (from 3c). Only residues named IPA/ACN/IMD/ACEA/PHOL/ACOH will be added.")
# ------------------------------------------------------------------------------

args = parser.parse_args()
suffix = "_no_ligand" if args.no_ligand else ""

# 1. Load receptor
receptor_pdb = PDBFile(args.input_receptor)

# 2. Force field setup
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 forcefield")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all forcefield")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# 3. Ligand logic (unchanged)
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

# ---- NEW: add Packmol probes if provided ------------------------------------
present_probe_resnames = set()
if args.packmol_probes_pdb is not None:
    print(f"🧪 [Packmol] Importing probes from: {args.packmol_probes_pdb}")
    pack_pdb = PDBFile(args.packmol_probes_pdb)
    probe_names = set(PROBE_DEF.keys())
    # Extract only probe residues
    sub = Modeller(pack_pdb.topology, pack_pdb.positions)
    # delete everything that is NOT one of the probe resnames
    to_delete = []
    for res in sub.topology.residues():
        if res.name not in probe_names:
            to_delete.extend(list(res.atoms()))
        else:
            present_probe_resnames.add(res.name)
    if len(to_delete) > 0:
        sub.delete(to_delete)
    if sum(1 for _ in sub.topology.atoms()) == 0:
        print("⚠️ Packmol file contains no residues named IPA/ACN/IMD/ACEA/PHOL/ACOH. Skipping.")
    else:
        # Register SMIRNOFF templates for molecules we might add
        off_mols = []
        for resname in sorted(present_probe_resnames):
            smiles = PROBE_DEF[resname]
            off = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
            off.generate_conformers(n_conformers=1)
            off_mols.append(off)
        if off_mols:
            smirnoff_all = SMIRNOFFTemplateGenerator(molecules=off_mols)
            forcefield.registerTemplateGenerator(smirnoff_all.generator)

        # Add all probes at their positions
        modeller.add(sub.topology, sub.positions)
        print(f"✅ Added probes: {', '.join(sorted(present_probe_resnames))}")
# -----------------------------------------------------------------------------

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
