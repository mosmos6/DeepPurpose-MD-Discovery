# 5_md_simulation.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Run OpenMM-based 2 ns MD simulation of ligand-receptor complex using OpenFF

from openmm.app import *
from openmm import MonteCarloBarostat, Platform, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm
from sys import stdout
import mdtraj as md
import matplotlib.pyplot as plt

# === Load Receptor ===
receptor_pdb = PDBFile("receptor_cleaned.pdb")

# === Load Ligand ===
ligand = Molecule.from_file("ligand.sdf")
ligand_positions = to_openmm(ligand.conformers[0])
ligand_top = OFFTopology.from_molecules([ligand]).to_openmm()

# === Define Force Field ===
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
smirnoff = SMIRNOFFTemplateGenerator(molecules=[ligand])
forcefield.registerTemplateGenerator(smirnoff.generator)

# === Merge Systems ===
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)
modeller.add(ligand_top, ligand_positions)
modeller.addHydrogens(forcefield)

with open("combined_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("✅ Receptor + ligand merged and hydrogens added.")

# === Solvation ===
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(7.0, 7.0, 7.0) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)

with open("solvated_receptor_ligand.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("✅ Solvated system ready.")

# === System Definition ===
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1 * nanometer,
    constraints=HBonds
)
system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin, 25))

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
platform = Platform.getPlatformByName("CUDA")
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# === Energy Minimization ===
print("🔹 Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter("minimized.pdb", 100))

# === NVT Equilibration ===
print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * kelvin)
simulation.reporters.append(PDBReporter("nvt_equilibrated.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

# === NPT Equilibration ===
print("🔹 NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter("npt_equilibrated.pdb", 500))
simulation.step(2500)  # 5 ps

# === Production MD (2 ns) ===
print("🔥 Production MD: 2 ns (1,000,000 steps)")
simulation.reporters.append(DCDReporter("production_md.dcd", 1000))
simulation.reporters.append(StateDataReporter("production_md.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(1000000)  # 2 ns

# === Save Final Frame ===
with open("final_structure.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
print("🎉 MD complete → final_structure.pdb saved.")


