# Author: Iori Mochizuki + patch
# Updated: 2025-10-22
# Description: Run docking with AutoDock Vina for one or many ligands.

import argparse, subprocess, re
from pathlib import Path

def _list_ligand_pdbqt():
    """Return list of PDBQT ligand paths in preferred order."""
    files = []
    # prefer explicit numbered ligands if present
    files += sorted(Path(".").glob("ligand_*.pdbqt"),
                    key=lambda p: (int(re.sub(r"[^\d]", "", p.stem) or 0), p.stem))
    # fall back to legacy single name
    if Path("ligand.pdbqt").exists():
        files.append(Path("ligand.pdbqt"))
    # de‑duplicate while preserving order
    seen, out = set(), []
    for p in files:
        if p.resolve() not in seen:
            out.append(p)
            seen.add(p.resolve())
    return out

def _compute_center(receptor_path: Path, args):
    # Priority 1: manual
    if args.center_x is not None and args.center_y is not None and args.center_z is not None:
        x = round(float(args.center_x), 3); y = round(float(args.center_y), 3); z = round(float(args.center_z), 3)
        print(f"📍 Docking box center (manual): x={x} y={y} z={z}")
        return x, y, z

    # Priority 2: RNA range
    def _parse_range(s):
        for dash in ['–', '—', '−', '‒', '―']: s = s.replace(dash, '-')
        a, b = s.split('-'); return int(a), int(b)

    if args.rna and args.res_range:
        start, end = _parse_range(args.res_range)
        s = [0.0, 0.0, 0.0]; n = 0
        for ln in receptor_path.read_text().splitlines():
            if ln.startswith("ATOM"):
                resid = int(ln[22:26])
                if start <= resid <= end:
                    s[0] += float(ln[30:38]); s[1] += float(ln[38:46]); s[2] += float(ln[46:54]); n += 1
        if n > 0:
            x, y, z = (round(s[0]/n,3), round(s[1]/n,3), round(s[2]/n,3))
            print(f"📍 RNA centroid {start}-{end}: x={x} y={y} z={z}")
            return x, y, z
        print("⚠️ No atoms in specified RNA range; falling back.")

    # Priority 3: HETATM centroid
    if args.use_residue_centroid:
        resnames = set()
        for ln in receptor_path.read_text().splitlines():
            if ln.startswith("HETATM"):
                resnames.add(ln[17:20].strip())
        s = [0.0,0.0,0.0]; n = 0
        for ln in receptor_path.read_text().splitlines():
            if ln.startswith("HETATM") and (ln[17:20].strip() in resnames):
                s[0]+=float(ln[30:38]); s[1]+=float(ln[38:46]); s[2]+=float(ln[46:54]); n+=1
        if n>0:
            x,y,z = (round(s[0]/n,3), round(s[1]/n,3), round(s[2]/n,3))
            print(f"📍 HETATM centroid: x={x} y={y} z={z}")
            return x,y,z
        print("⚠️ No HETATM residues; falling back.")

    # Priority 4: all ATOM centroid
    s = [0.0,0.0,0.0]; n = 0
    for ln in receptor_path.read_text().splitlines():
        if ln.startswith("ATOM"):
            s[0]+=float(ln[30:38]); s[1]+=float(ln[38:46]); s[2]+=float(ln[46:54]); n += 1
    x, y, z = (round(s[0]/n,3), round(s[1]/n,3), round(s[2]/n,3))
    print(f"📍 ATOM centroid: x={x} y={y} z={z}")
    return x, y, z

def _extract_model1_pdb(in_pdbqt: Path, out_pdb: Path):
    inside = False
    with in_pdbqt.open("r") as fin, out_pdb.open("w") as fout:
        for line in fin:
            if line.startswith("MODEL 2"): break
            if line.startswith("MODEL 1"): inside = True
            if inside and not line.startswith("MODEL"):
                # keep only (HET)ATOM,TER lines (let Step 4 clean elements)
                if line.startswith(("ATOM","HETATM","TER")):
                    fout.write(line)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--use-residue-centroid", action="store_true")
    p.add_argument("--rna", action="store_true")
    p.add_argument("--res-range", type=str, default=None)
    p.add_argument("--center_x", type=float, default=None)
    p.add_argument("--center_y", type=float, default=None)
    p.add_argument("--center_z", type=float, default=None)
    p.add_argument("--vina_seed", type=int, default=12345)
    args = p.parse_args()

    receptor_for_center = Path("receptor_for_centroid.pdb")
    receptor_pdbqt = Path("receptor.pdbqt")
    vina_bin = "./vina_1.2.5_linux_x86_64"

    ligands = _list_ligand_pdbqt()
    if not ligands:
        raise SystemExit("No ligand PDBQT files found (looked for ligand_*.pdbqt or ligand.pdbqt).")

    cx, cy, cz = _compute_center(receptor_for_center, args)

    print(f"Found {len(ligands)} ligand(s) → {', '.join(p.name for p in ligands)}")

    # If only one ligand and named exactly 'ligand.pdbqt', keep legacy filenames too
    legacy_single = (len(ligands) == 1 and ligands[0].name == "ligand.pdbqt")

    for i, lig in enumerate(ligands, start=1):
        tag = lig.stem             # 'ligand' or 'ligand_3'
        out_pdbqt = Path(f"output_{tag}.pdbqt")
        out_best  = Path(f"best_docked_{tag}.pdb")
        out_clean = Path(f"output_{tag}.pdb")

        cmd = (
            f"{vina_bin} "
            f"--receptor {receptor_pdbqt} "
            f"--ligand {lig} "
            f"--out {out_pdbqt} "
            f"--center_x {cx} --center_y {cy} --center_z {cz} "
            f"--size_x 30 --size_y 30 --size_z 30 "
            f"--seed {int(args.vina_seed)} --exhaustiveness 20"
        )
        print(f"\n▶️  [{i}/{len(ligands)}] Vina docking for {lig.name}")
        subprocess.run(cmd, shell=True, check=True)
        _extract_model1_pdb(out_pdbqt, out_best)
        # final cleaned pose used by Step 4
        out_clean.write_text(out_best.read_text())
        print(f"✅ {tag}: best pose → {out_best}  | cleaned → {out_clean}")

        if legacy_single:
            # also write legacy filenames expected by old scripts
            Path("output.pdbqt").write_bytes(out_pdbqt.read_bytes())
            Path("best_docked_ligand.pdb").write_bytes(out_best.read_bytes())
            Path("output.pdb").write_bytes(out_clean.read_bytes())
            print("↪️  Also wrote legacy output.pdbqt / best_docked_ligand.pdb / output.pdb")

if __name__ == "__main__":
    main()
