# 4_align_ligand.py
# Author: Iori Mochizuki + minimal multi-ligand loop
# Updated: 2025-10-30
# Description:
#   Iterate the original single‑ligand alignment for all names in ligands.index.
#   Steps per ligand: Kabsch align → write aligned PDB → fix element column →
#   convert to SDF via Open Babel. Produces ligands_sdf.index for Step 5.

from pathlib import Path
import numpy as np
import subprocess

from utils import extract_coordinates, kabsch, fix_pdb_element_column

INDEX_FILE = Path("ligands.index")

def _read_names() -> list[str]:
    if not INDEX_FILE.exists():
        raise SystemExit("ligands.index not found (run Step 1).")
    names = []
    seen = set()
    for ln in INDEX_FILE.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if s not in seen:
            names.append(s)
            seen.add(s)
    if not names:
        raise SystemExit("ligands.index is empty.")
    return names

def _align_one(name: str):
    # --- file layout identical to your original single‑ligand step ---
    docked_pdb      = Path(f"output_{name}.pdb")
    original_pdb    = Path(f"{name}.pdb")
    stripped_pdbqt  = Path(f"{name}.pdbqt")
    aligned_pdb     = Path(f"aligned_{name}_fixed.pdb")
    final_pdb       = Path(f"fixed_{name}.pdb")
    final_sdf       = Path(f"{name}.sdf")

    # minimal sanity checks
    missing = [p.name for p in (docked_pdb, original_pdb, stripped_pdbqt) if not p.exists()]
    if missing:
        print(f"⚠️  Skipping {name}: missing {', '.join(missing)}")
        return False

    # === Step 1: Load coordinates (same helpers as before) ===
    docked_coords, _        = extract_coordinates(str(docked_pdb))
    full_coords,  full_atoms = extract_coordinates(str(original_pdb))
    stripped_coords, _      = extract_coordinates(str(stripped_pdbqt))

    print(f"✅ Loaded {len(docked_coords)} docked, {len(full_coords)} full, {len(stripped_coords)} stripped atoms for {name}")

    # === Step 2: Sort by radial distance (stability) ===
    d1 = np.linalg.norm(docked_coords   - docked_coords.mean(axis=0), axis=1)
    d2 = np.linalg.norm(stripped_coords - stripped_coords.mean(axis=0), axis=1)
    sorted_docked   = docked_coords[np.argsort(d1)]
    sorted_stripped = stripped_coords[np.argsort(d2)]

    # === Step 3: Kabsch alignment ===
    R, c_stripped, c_docked = kabsch(sorted_stripped, sorted_docked)
    aligned_coords = (full_coords - c_stripped) @ R.T + c_docked  # Å

    # === Step 4: Write aligned structure ===
    with aligned_pdb.open("w") as f:
        for i, (_, _, _, _, line) in enumerate(full_atoms):
            x, y, z = aligned_coords[i]
            f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")

    print(f"✅ Aligned ligand written to: {aligned_pdb}")

    # === Step 5: Fix element column ===
    fix_pdb_element_column(str(aligned_pdb), str(final_pdb))

    # === Step 6: Convert to SDF via Open Babel (unchanged route) ===
    try:
        # Pass Path objects directly, just like your original flow allowed
        subprocess.run(["obabel", final_pdb, "-O", final_sdf], check=True)
    except FileNotFoundError:
        raise SystemExit("Open Babel ('obabel') not found in PATH. Install it in your conda env.")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Open Babel failed converting {final_pdb} → {final_sdf}.\n{e}")

    print(f"✅ {name}: aligned → {aligned_pdb} | fixed → {final_pdb} | SDF → {final_sdf}")
    return True

def main():
    names = _read_names()
    print(f"\n[4] Aligning {len(names)} ligand(s): {', '.join(names)}\n")

    done = []
    for n in names:
        if _align_one(n):
            done.append(n)

    if not done:
        raise SystemExit("No ligands processed. Check inputs listed in ligands.index.")

    # Write the SDF index **as stems** (no .sdf), for Step 5 consumption
    Path("ligands_sdf.index").write_text("\n".join(done) + "\n", encoding="utf-8")

    print("\nSummary:")
    print(" • ligands_sdf.index (one name per line; SDF is <name>.sdf)")
    for n in done:
        print(f"   - {n}.sdf")

if __name__ == "__main__":
    main()
