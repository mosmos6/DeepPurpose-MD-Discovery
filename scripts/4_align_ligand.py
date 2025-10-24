# Author: Iori Mochizuki + patch
# Updated: 2025-10-22
# Description: Align one or many ligands to docked coordinates; export SDF for OpenFF.

import re, subprocess
from pathlib import Path
import numpy as np
from utils import extract_coordinates, kabsch, fix_pdb_element_column

def _discover_pairs():
    """
    Return list of (base_tag, files_dict) where base_tag is 'ligand' or 'ligand_3',
    and files_dict contains keys: original_pdb, stripped_pdbqt, docked_pdb.
    """
    tags = set()
    # prefer numbered
    for p in Path(".").glob("ligand_*.pdb"):
        tags.add(p.stem)
    for p in Path(".").glob("ligand_*.pdbqt"):
        tags.add(p.stem)
    for p in Path(".").glob("output_ligand_*.pdb"):
        m = re.match(r"output_(ligand_\d+)$", p.stem)
        if m: tags.add(m.group(1))
    # legacy single
    if Path("ligand.pdb").exists() or Path("ligand.pdbqt").exists() or Path("output.pdb").exists():
        tags.add("ligand")

    out = []
    for tag in sorted(tags, key=lambda t: (0 if t=="ligand" else 1, int(re.sub(r"[^\d]","",t) or 0))):
        files = {
            "original_pdb": Path(f"{tag}.pdb") if Path(f"{tag}.pdb").exists() else Path("ligand.pdb"),
            "stripped_pdbqt": Path(f"{tag}.pdbqt") if Path(f"{tag}.pdbqt").exists() else Path("ligand.pdbqt"),
            "docked_pdb": Path(f"output_{tag}.pdb") if Path(f"output_{tag}.pdb").exists() else Path("output.pdb"),
        }
        if files["original_pdb"].exists() and files["stripped_pdbqt"].exists() and files["docked_pdb"].exists():
            out.append((tag, files))
    return out

def _align_one(tag, files):
    docked_coords, docked_atoms = extract_coordinates(str(files["docked_pdb"]))
    full_coords,   full_atoms   = extract_coordinates(str(files["original_pdb"]))
    stripped_coords,_           = extract_coordinates(str(files["stripped_pdbqt"]))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {tag}")

    # sort by radial distance to centroid for stability (as in your original)
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned_coords = (full_coords - c_stripped) @ R.T + c_docked

    aligned_pdb = Path(f"aligned_{tag}_fixed.pdb")
    final_pdb   = Path(f"fixed_{tag}.pdb")
    final_sdf   = Path(f"{tag}.sdf")

    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned_coords[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    fix_pdb_element_column(str(aligned_pdb), str(final_pdb))
    subprocess.run(["obabel", str(final_pdb), "-O", str(final_sdf)], check=True)
    print(f"✅ {tag}: aligned → {aligned_pdb} | fixed → {final_pdb} | SDF → {final_sdf}")

    # legacy fallbacks if the tag is exactly 'ligand'
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
