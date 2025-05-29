# 5b_md_analysis_RMSD_RMSF.py
# Author: Iori Mochizuki
# Description: Analyze RMSD and RMSF after MD simulation

import mdtraj as md
import matplotlib.pyplot as plt

# === Load trajectory ===
print("📥 Loading trajectory and reference...")
traj = md.load("production_md.dcd", top="solvated_receptor_ligand.pdb")
ref = md.load("npt_equilibrated.pdb")

# === Analysis: RMSD ===
print("📊 RMSD analysis...")
rmsd = md.rmsd(traj, ref, frame=0)

plt.figure()
plt.plot(rmsd)
plt.xlabel("Frame")
plt.ylabel("RMSD (nm)")
plt.title("Ligand RMSD Over Time")
plt.tight_layout()
plt.savefig("quick_rmsd_plot.png")
print("✅ RMSD plot saved as quick_rmsd_plot.png")

# === Analysis: RMSF ===
print("📊 RMSF analysis...")
ca_atoms = traj.topology.select("name CA")
rmsf = md.rmsf(traj, ref, atom_indices=ca_atoms)

plt.figure()
plt.plot(rmsf)
plt.xlabel("Residue Index")
plt.ylabel("RMSF (nm)")
plt.title("Cα RMSF Profile")
plt.tight_layout()
plt.savefig("quick_rmsf_plot.png")
print("✅ RMSF plot saved as quick_rmsf_plot.png")
