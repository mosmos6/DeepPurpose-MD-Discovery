# scripts/4_align_ligand.py
# Author: Iori Mochizuki + patch
# Updated: 2025-10-29
# Description: Align one or many ligands to docked coordinates; export SDF for OpenFF.
# Key change: write SDF with RDKit from the original Step-1 MOL2 topology
#             to preserve bond orders & formal charges (no Babel ambiguity).

import re
from pathlib import Path
import numpy as np
from utils import extract_coordinates, kabsch, fix_pdb_element_column

# NEW: RDKit is used here to write SDF using original MOL2 topology
from rdkit import Chem

def _discover_pairs():
    """
    Return list of (tag, files) where tag is 'ligand' or 'ligand_3',
    files has: original_pdb, stripped_pdbqt, docked_pdb, template_mol2.
    """
    tags = set()
    for p in Path(".").glob("ligand_*.pdb"):
        tags.add(p.stem)
    for p in Path(".").glob("ligand_*.pdbqt"):
        tags.add(p.stem)
    for p in Path(".").glob("output_ligand_*.pdb"):
        m = re.match(r"output_(ligand_\d+)$", p.stem)
        if m: tags.add(m.group(1))
    if Path("ligand.pdb").exists() or Path("ligand.pdbqt").exists() or Path("output.pdb").exists():
        tags.add("ligand")

    def _order_key(t):
        return (0 if t == "ligand" else 1, int(re.sub(r"[^\d]", "", t) or 0))

    out = []
    for tag in sorted(tags, key=_order_key):
        files = {
            "original_pdb": Path(f"{tag}.pdb") if Path(f"{tag}.pdb").exists() else Path("ligand.pdb"),
            "stripped_pdbqt": Path(f"{tag}.pdbqt") if Path(f"{tag}.pdbqt").exists() else Path("ligand.pdbqt"),
            "docked_pdb": Path(f"output_{tag}.pdb") if Path(f"output_{tag}.pdb").exists() else Path("output.pdb"),
            # NEW: template MOL2 from Step 1 (same atom ordering as original PDB)
            "template_mol2": Path(f"{tag}.mol2") if Path(f"{tag}.mol2").exists() else Path("ligand.mol2"),
        }
        if all(files[k].exists() for k in ("original_pdb", "stripped_pdbqt", "docked_pdb")):
            out.append((tag, files))
    return out

def _write_sdf_from_mol2_with_coords(tag: str, mol2_path: Path, coords_xyz: np.ndarray, out_sdf: Path):
    """
    Load RDKit Mol from MOL2 (preserves bond orders & charges),
    replace conformer coordinates with aligned coords, and write SDF.
    """
    # Load MOL2 as RDKit mol (keeps charges/bonds)
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    if mol is None:
        raise RuntimeError(f"[{tag}] RDKit failed to read MOL2: {mol2_path}")

    n = mol.GetNumAtoms()
    if n != coords_xyz.shape[0]:
        raise RuntimeError(f"[{tag}] Atom count mismatch: MOL2 has {n}, aligned coords have {coords_xyz.shape[0]}")

    # Ensure a conformer exists; create if missing
    if mol.GetNumConformers() == 0:
        conf = Chem.Conformer(n)
        mol.AddConformer(conf, assignId=True)
    conf = mol.GetConformer(0)

    # Set coordinates (Å)
    for i in range(n):
        x, y, z = float(coords_xyz[i,0]), float(coords_xyz[i,1]), float(coords_xyz[i,2])
        conf.SetAtomPosition(i, (x, y, z))

    # Write SDF
    writer = Chem.SDWriter(str(out_sdf))
    writer.write(mol)
    writer.close()

def _align_one(tag, files):
    docked_coords, docked_atoms = extract_coordinates(str(files["docked_pdb"]))
    full_coords,   full_atoms   = extract_coordinates(str(files["original_pdb"]))
    stripped_coords,_           = extract_coordinates(str(files["stripped_pdbqt"]))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {tag}")

    # sort by radial distance to centroid for stability
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned_coords = (full_coords - c_stripped) @ R.T + c_docked  # Å

    aligned_pdb = Path(f"aligned_{tag}_fixed.pdb")
    final_pdb   = Path(f"fixed_{tag}.pdb")
    final_sdf   = Path(f"{tag}.sdf")

    # Write aligned PDB (for manual QC) and fix element column
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned_coords[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
    fix_pdb_element_column(str(aligned_pdb), str(final_pdb))

    # NEW: Write SDF from original MOL2 topology but with aligned coordinates
    if not files["template_mol2"].exists():
        raise FileNotFoundError(f"[{tag}] Missing template MOL2 from Step 1: {files['template_mol2']}")
    _write_sdf_from_mol2_with_coords(tag, files["template_mol2"], aligned_coords, final_sdf)

    print(f"✅ {tag}: aligned → {aligned_pdb} | fixed → {final_pdb} | SDF (RDKit) → {final_sdf}")

    # legacy single-ligand copies, if needed
    if tag == "ligand":
        Path("aligned_ligand_fixed.pdb").write_bytes(aligned_pdb.read_bytes())
        Path("fixed_ligand.pdb").write_bytes(final_pdb.read_bytes())
        Path("ligand.sdf").write_bytes(final_sdf.read_bytes())
        print("↪️  Also wrote legacy aligned_ligand_fixed.pdb / fixed_ligand.pdb / ligand.sdf")

def main():
    pairs = _discover_pairs()
    if not pairs:
        raise SystemExit("No ligand triplets found (expected ligand[_N].pdb/.pdbqt and output_ligand[_N].pdb).")
    for tag, files in pairs:
        print(f"\n--- Aligning {tag} ---")
        _align_one(tag, files)

if __name__ == "__main__":
    main()
