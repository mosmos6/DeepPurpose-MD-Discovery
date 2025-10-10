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

# --- MixMD-lite utils (kept separate to avoid math/unit collisions) ---
import os, sys
sys.path.append(os.path.dirname(__file__))  # allow "scripts/..." imports when run from repo root
from mixmd_utils_5c import (
    PROBES,
    parse_probe_fractions,
    build_probe_molecules,
    register_smirnoff_templates,
    replace_waters_with_probes,
)

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

# 1. Load receptor
receptor_pdb = PDBFile(args.input_receptor)

# 2. Force field setup
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 forcefield")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all forcefield")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# --- Probe setup (only when --mixmd is on) ---
probe_fracs = parse_probe_fractions(args.probe_fractions) if args.mixmd else {}
unique_off_mols = []

# include the “main ligand” only if not in apo mode
if not args.no_ligand:
    ligand = Molecule.from_file(args.input_ligand)
    unique_off_mols.append(ligand)

probe_cache = {}
if args.mixmd:
    # build OFF molecules for all requested probes (1 conformer each, centered)
    probe_cache = build_probe_molecules(list(probe_fracs.keys()))
    # register a single SMIRNOFF template for ligand(+probes) before creating the system
    unique_off_mols.extend([v["off"] for v in probe_cache.values()])

if unique_off_mols:
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

if args.mixmd:
    print("🧪 [MixMD-lite] Replacing waters with probe molecules...")
    replace_waters_with_probes(
        modeller,
        probe_cache=probe_cache,
        probe_fracs=probe_fracs,
        total_count=int(args.probe_total),
        seed=int(args.probe_seed),
    )
    # (optional) write a snapshot after placement so you can visually inspect
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
