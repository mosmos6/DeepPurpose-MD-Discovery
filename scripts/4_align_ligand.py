# scripts/4_align_ligand.py
# Author: Iori Mochizuki + patch
# Updated: 2025-10-28
# Description: Align one or many ligands to their docked poses; export SDF for OpenFF.
# Inputs expected (per ligand tag <T>):
#   <T>.pdb       (original RDKit build from Step 1)
#   <T>.pdbqt     (stripped used for docking)
#   output_<T>.pdb  (MODEL 1 best pose extracted by Step 3)
#
# Outputs (per <T>):
#   aligned_<T>_fixed.pdb, fixed_<T>.pdb, <T>.sdf
# Plus:
#   ligands_sdf.index  (list of produced SDFs, one per line)
#   Legacy copies if there is exactly one ligand and its tag != 'ligand'

from pathlib import Path
import re, subprocess, sys
import numpy as np
from utils import extract_coordinates, kabsch, fix_pdb_element_column

def _load_ligand_tags():
    """
    Preferred: read ligands.index (one name per line) written by Step 1.
    Fallback: discover by scanning *.pdbqt (use basenames).
    """
    idx = Path("ligands.index")
    tags = []
    if idx.exists():
        for ln in idx.read_text(encoding="utf-8").splitlines():
            t = ln.strip()
            if t:
                tags.append(t)
    else:
        # fallback discovery
        for p in sorted(Path(".").glob("*.pdbqt")):
            t = p.stem
            if t not in tags:
                tags.append(t)
    # keep original order; de-dup
    seen, out = set(), []
    for t in tags:
        if t not in seen:
            out.append(t); seen.add(t)
    return out

def _files_for_tag(tag: str):
    files = {
        "original_pdb": Path(f"{tag}.pdb"),
        "stripped_pdbqt": Path(f"{tag}.pdbqt"),
        "docked_pdb": Path(f"output_{tag}.pdb"),
    }
    missing = [k for k, p in files.items() if not p.exists()]
    return files, missing

def _align_one(tag: str, files: dict):
    docked_coords, docked_atoms = extract_coordinates(str(files["docked_pdb"]))
    full_coords,   full_atoms   = extract_coordinates(str(files["original_pdb"]))
    stripped_coords,_           = extract_coordinates(str(files["stripped_pdbqt"]))

    print(f"✅ {tag}: loaded  docked={len(docked_coords)}  full={len(full_coords)}  stripped={len(stripped_coords)}")

    # Sort by radial distance to centroids for numerical stability (as in your original)
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    # Kabsch on stripped→docked; apply to full
    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned_coords = (full_coords - c_stripped) @ R.T + c_docked

    aligned_pdb = Path(f"aligned_{tag}_fixed.pdb")
    final_pdb   = Path(f"fixed_{tag}.pdb")
    final_sdf   = Path(f"{tag}.sdf")

    # Write aligned PDB (preserve original record formatting)
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned_coords[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    fix_pdb_element_column(str(aligned_pdb), str(final_pdb))

    # PDB → SDF (Open Babel)
    subprocess.run(["obabel", str(final_pdb), "-O", str(final_sdf)], check=True)
    print(f"   ↳ aligned → {aligned_pdb} | fixed → {final_pdb} | SDF → {final_sdf}")

    return final_sdf

def main():
    tags = _load_ligand_tags()
    if not tags:
        sys.exit("No ligands found. Expecting ligands.index or *.pdbqt from Step 1.")

    produced_sdfs = []
    print(f"\n[4] Aligning {len(tags)} ligand(s): {', '.join(tags)}\n")

    for tag in tags:
        files, missing = _files_for_tag(tag)
        if missing:
            print(f"⚠️  Skipping {tag} (missing: {', '.join(missing)})")
            continue
        sdf = _align_one(tag, files)
        produced_sdfs.append(sdf)

    if not produced_sdfs:
        sys.exit("No ligands were aligned (see warnings above).")

    # Index of produced SDFs for MD
    Path("ligands_sdf.index").write_text(
        "\n".join(str(p) for p in produced_sdfs) + "\n",
        encoding="utf-8"
    )
    print("\n✅ Wrote ligands_sdf.index")

    # Legacy compatibility when exactly one ligand and tag != 'ligand'
    if len(tags) == 1 and tags[0] != "ligand":
        tag = tags[0]
        # copy outputs to legacy names
        Path("aligned_ligand_fixed.pdb").write_bytes(Path(f"aligned_{tag}_fixed.pdb").read_bytes())
        Path("fixed_ligand.pdb").write_bytes(Path(f"fixed_{tag}.pdb").read_bytes())
        Path("ligand.sdf").write_bytes(Path(f"{tag}.sdf").read_bytes())
        print("↪️  Also wrote legacy aligned_ligand_fixed.pdb / fixed_ligand.pdb / ligand.sdf")

if __name__ == "__main__":
    main()
