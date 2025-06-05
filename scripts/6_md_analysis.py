# 6_md_analysis.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Description: Analyze MD trajectory (ligand binding, RMSD/F, PCA, FEL)

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.decomposition import PCA
from scipy.constants import Boltzmann

# === Input paths ===
traj_path = "production_md.dcd"
top_path = "final_structure.pdb"
ligand_resname = "UNK"  # Update if your ligand uses a different residue name

# === Load trajectory ===
traj = md.load_dcd(traj_path, top=top_path)

# -------------------------------
# 1. Ligand RMSD Over Time
# -------------------------------
ligand_atoms = traj.topology.select(f"resname {ligand_resname}")
ligand_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=ligand_atoms)

plt.figure(figsize=(8, 5))
plt.plot(traj.time / 1000, ligand_rmsd)
plt.xlabel("Time (ns)")
plt.ylabel("Ligand RMSD (nm)")
plt.title("Ligand RMSD Over Time")
plt.grid()
plt.savefig("Ligand_RMSD.png")
plt.close()
print("✅ Ligand RMSD saved.")

# -------------------------------
# 2. Ligand-Residue Interaction Frequency (%)
# -------------------------------
protein_atoms = traj.topology.select("protein")
atom_pairs = np.array([(lig, prot) for lig in ligand_atoms for prot in protein_atoms])
distances = md.compute_distances(traj, atom_pairs)
threshold = 0.4  # 4 Å

contact_frames = distances < threshold
residue_contact_counts = defaultdict(int)

for frame in contact_frames:
    for idx, contact in enumerate(frame):
        if contact:
            res = traj.topology.atom(atom_pairs[idx][1]).residue.name
            residue_contact_counts[res] += 1

total_contacts = sum(residue_contact_counts.values())
freq_dict = {res: (cnt / total_contacts) * 100 for res, cnt in residue_contact_counts.items()}
freq_sorted = sorted(freq_dict.items(), key=lambda x: x[1], reverse=True)

print("🔹 Residue-Ligand Interaction Frequencies (%):")
for res, freq in freq_sorted:
    print(f"{res}: {freq:.2f}%")

# -------------------------------
# 3. Most Frequent Residue Contacts
# -------------------------------
lig_resid = {a.residue.index for a in traj.topology.atoms if a.residue.name == ligand_resname}
prot_resids = {a.residue.index for a in traj.topology.atoms if a.residue.is_protein}
pairs = np.array([(l, p) for l in lig_resid for p in prot_resids])
contacts, _ = md.compute_contacts(traj, pairs)

counts = np.sum(contacts < threshold, axis=0)
res_freq = {}
for i, c in enumerate(counts):
    if c > 0:
        res = traj.topology.residue(pairs[i, 1])
        res_freq[res] = c

for res, count in sorted(res_freq.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"{res}: {count} frames")

# -------------------------------
# 4. Ligand Residence Time near Target Residue
# -------------------------------
target_resid = 124  # e.g. LYS125 → MDTraj uses 0-indexed
target_atoms = traj.topology.select(f"resname LYS and resid {target_resid}")

if target_atoms.any():
    ap = np.array([(l, t) for l in ligand_atoms for t in target_atoms])
    d = md.compute_distances(traj, ap)
    bound = np.any(d < threshold, axis=1)
    percent = (np.sum(bound) / traj.n_frames) * 100
    print(f"✅ Residence near LYS125: {percent:.2f}% of trajectory")
else:
    print("⚠️ LYS125 not found. Adjust residue index.")

# -------------------------------
# 5. Hydration Shell Around Ligand
# -------------------------------
water_atoms = traj.topology.select("water")
water_pairs = np.array([(l, w) for l in ligand_atoms for w in water_atoms])
wd = md.compute_distances(traj, water_pairs)
hydration_counts = np.sum(wd < 0.35, axis=1)

plt.figure()
plt.plot(traj.time / 1000, hydration_counts)
plt.xlabel("Time (ns)")
plt.ylabel("Water Contacts (<3.5Å)")
plt.title("Ligand Hydration Over Time")
plt.savefig("Ligand_Hydration.png")
plt.close()
print("✅ Hydration plot saved.")

# -------------------------------
# 6. Principal Component Analysis (PCA)
# -------------------------------
traj.superpose(traj, 0)
backbone = traj.topology.select("name CA")
xyz_flat = traj.xyz[:, backbone, :].reshape(traj.n_frames, -1)

pca = PCA(n_components=2)
pcs = pca.fit_transform(xyz_flat)

plt.figure()
plt.plot(pcs[:, 0], pcs[:, 1], alpha=0.7)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("Protein PCA Motion")
plt.grid(True)
plt.savefig("PCA_Motion.png")
plt.close()
print("✅ PCA projection saved.")

# -------------------------------
# 7. Free Energy Landscape (FEL)
# -------------------------------
def compute_fel(pc1, pc2, bins=100):
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins, density=True)
    P = H.T
    P[P == 0] = 1e-12
    T = 310.15  # body temp K
    G = -Boltzmann * T * np.log(P)
    return G, xedges, yedges

def plot_fel(G, xedges, yedges, title, filename):
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    plt.figure(figsize=(7, 5))
    plt.contourf(X, Y, G, levels=30, cmap="viridis")
    plt.colorbar(label="Free Energy (J)")
    plt.title(title)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

G, xe, ye = compute_fel(pcs[:, 0], pcs[:, 1])
plot_fel(G, xe, ye, "Free Energy Landscape (FEL)", "FEL_Ligand.png")
print("✅ FEL plot saved.")
