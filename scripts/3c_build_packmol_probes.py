# scripts/3c_build_packmol_probes.py
# Universal Packmol driver for MixMD-style placement (protein/RNA targets, incl. 40S)
# Author: Iori Mochizuki (pipeline owner) + collaborator patch
# Updated: 2025-10-14

import argparse, os, shutil, subprocess, sys, math, csv
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass

# ---- minimal RDKit-only helpers (no OpenMM/OpenFF here)
from rdkit import Chem
from rdkit.Chem import AllChem

# --- Box helper --------------------------------------------------------------
from dataclasses import dataclass

@dataclass(frozen=True)
class Box:
    """Simple cubic box (nm)."""
    side_nm: float

    @property
    def half_nm(self) -> float:
        return 0.5 * self.side_nm

    @property
    def lo_nm(self) -> float:
        return -self.half_nm

    @property
    def hi_nm(self) -> float:
        return self.half_nm
# ---------------------------------------------------------------------------

PROBES: Dict[str, Tuple[str, str]] = {
    # key: (RESNAME, SMILES)
    "ipa":   ("IPA",  "CC(C)O"),        # isopropanol
    "acn":   ("ACN",  "CC#N"),          # acetonitrile
    "imd":   ("IMD",  "c1[nH]cnc1"),    # imidazole (neutral, N1 protonated)
    "aceam": ("ACEA", "CC(=O)N"),       # acetamide
    "phol":  ("PHOL", "c1ccc(cc1)O"),   # phenol
    "acoh":  ("ACOH", "CC(=O)O"),       # acetic acid
}

@dataclass
class Box:
    L: float  # nm
    lo: float # nm
    hi: float # nm

def parse_kv_fracs(s: str) -> Dict[str, float]:
    fr = {}
    for tok in s.split(","):
        tok = tok.strip()
        if not tok: 
            continue
        if "=" not in tok:
            raise ValueError(f"Bad token in --probe-fractions: '{tok}' (use key=val)")
        k, v = tok.split("=", 1)
        k = k.strip().lower()
        if k not in PROBES:
            raise ValueError(f"Unknown probe key '{k}'. Allowed: {list(PROBES.keys())}")
        fr[k] = float(v.strip())
    if not fr:
        raise ValueError("At least one probe fraction must be provided.")
    ssum = sum(fr.values())
    if ssum <= 0:
        raise ValueError("Sum of fractions must be > 0.")
    for k in fr:
        fr[k] /= ssum
    return fr

def rdkit_centered_pdb(mol: Chem.Mol, resname: str) -> str:
    """Return a PDB block with the molecule centered at its centroid and residue name set."""
    # Ensure 3D conformer
    mol = Chem.AddHs(mol)
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    coords = []
    cx = cy = cz = 0.0
    for i in range(n):
        x, y, z = conf.GetAtomPosition(i)
        coords.append((x, y, z))  # Å
        cx += x; cy += y; cz += z
    cx /= n; cy /= n; cz /= n
    # Translate to centroid 0, then write PDB lines ourselves to control resname
    lines = []
    for idx, (x, y, z) in enumerate(coords, start=1):
        x0 = x - cx; y0 = y - cy; z0 = z - cz
        atom = mol.GetAtomWithIdx(idx-1)
        elem = atom.GetSymbol().rjust(2)
        name = (elem + "  ").ljust(4)[:4]  # crude but stable
        # HETATM serial, name, resname, chain, resid, x,y,z, occ, b, element
        lines.append(f"HETATM{idx:5d} {name} {resname:>3s} Z{1:4d}    "
                     f"{x0:8.3f}{y0:8.3f}{z0:8.3f}  1.00  0.00          {elem:>2s}")
    # TER + END
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"

def write_probe_templates(out_dir: Path) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for key, (resname, smi) in PROBES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise RuntimeError(f"SMILES parse failed for {key}: {smi}")
        pdb_block = rdkit_centered_pdb(mol, resname)
        path = out_dir / f"{resname}.pdb"
        path.write_text(pdb_block, encoding="utf-8")
        out[key] = path
    return out

def auto_box_from_receptor(receptor_pdb: Path, padding_nm: float) -> Box:
    """
    Infer a cubic box that encloses the receptor with margin = padding_nm (nm).
    PDB coords are in Å; we convert to nm internally.
    Returns: Box(side_nm=<float>)
    """
    coords_A = []
    with open(receptor_pdb, "r") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                # PDB columns: x[30:38], y[38:46], z[46:54] in Å
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords_A.append((x, y, z))
                except ValueError:
                    continue

    if not coords_A:
        raise ValueError("No ATOM/HETATM coordinates found in receptor PDB.")

    import numpy as _np
    xyz_nm = _np.asarray(coords_A, dtype=_np.float64) * 0.1  # Å → nm
    mins = xyz_nm.min(axis=0)
    maxs = xyz_nm.max(axis=0)
    span_nm = float((maxs - mins).max())
    side_nm = span_nm + 2.0 * float(padding_nm)

    # nice printout (unchanged from your logs)
    half = 0.5 * side_nm
    print(f"[3c] Box: {side_nm:.2f} nm (lo={-half:.2f}, hi={half:.2f})")

    return Box(side_nm=side_nm)

def write_packmol_input(inp_path: Path,
                        receptor_pdb: Path,
                        probe_pdbs: dict,
                        counts: dict,
                        box_side_nm: float,   # <— numeric (nm)
                        tol_A: float,
                        out_pdb: Path):
    """
    Emit a Packmol input file that:
      - fixes the receptor at the box center,
      - distributes probes uniformly inside the full box,
      - uses correct 'structure ... end structure' blocks.

    Units: Packmol expects Å. We convert nm -> Å where needed.
    """
    side_A = float(box_side_nm) * 10.0      # nm -> Å
    half_A = 0.5 * side_A
    lo, hi = -half_A, half_A

    lines = []
    # Global header
    lines.append(f"tolerance {float(tol_A):.3f}")  # Å
    lines.append("filetype pdb")
    lines.append("add_box_sides")
    lines.append(f"output {out_pdb.as_posix()}")

    # Receptor fixed at center
    lines.append(f"structure {receptor_pdb.as_posix()}")
    lines.append("  number 1")
    lines.append("  center")
    lines.append("  fixed 0. 0. 0. 0. 0. 0.")
    lines.append("end structure")

    # Probes across full orthorhombic box
    for key, n in counts.items():
        n_int = int(n)
        if n_int <= 0:
            continue
        probe_path = probe_pdbs[key]
        lines.append(f"structure {probe_path.as_posix()}")
        lines.append(f"  number {n_int}")
        lines.append(f"  inside box {lo:.3f} {lo:.3f} {lo:.3f} {hi:.3f} {hi:.3f} {hi:.3f}")
        lines.append("end structure")

    inp_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def run_packmol(input_file: Path):
    """
    Run Packmol by giving it a *seekable* stdin (the input file handle).
    This avoids the Fortran 'Illegal seek' when stdin is a pipe.

    Equivalent CLI:  packmol < input_file
    """
    exe = shutil.which("packmol")
    if exe is None:
        raise RuntimeError("Packmol executable not found in PATH.")

    print(f"[3c] Running: {exe} < {input_file}")

    # IMPORTANT: open a real file handle for stdin so Fortran REWIND works
    with open(input_file, "rb") as fin:
        proc = subprocess.run(
            [exe],
            stdin=fin,                      # <-- key change (no 'input=...')
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    if proc.returncode != 0:
        # show Packmol's own messages to help diagnose
        sys.stderr.write(proc.stdout.decode("utf-8", errors="ignore"))
        sys.stderr.write(proc.stderr.decode("utf-8", errors="ignore"))
        raise RuntimeError("[3c] Packmol failed. See console above.")

    print("[3c] Packmol completed.")


def parse_probe_centroids(packmol_pdb: Path, probe_resnames: List[str], out_csv: Path):
    # centroid per (resname, resid). PDB coords are Å.
    groups: Dict[Tuple[str,int], List[Tuple[float,float,float]]] = {}
    with open(packmol_pdb, "r") as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue
            resname = line[17:20].strip()
            if resname not in probe_resnames:
                continue
            resid = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])  # Å
            groups.setdefault((resname, resid), []).append((x, y, z))
    rows = []
    for (resname, resid), xyzs in groups.items():
        n = len(xyzs)
        cx = sum(p[0] for p in xyzs)/n * 0.1  # nm
        cy = sum(p[1] for p in xyzs)/n * 0.1
        cz = sum(p[2] for p in xyzs)/n * 0.1
        rows.append((resname, resid, cx, cy, cz))
    with open(out_csv, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["resname","resid","x_nm","y_nm","z_nm"])
        w.writerows(rows)
    print(f"[3c] Wrote centroids for {len(rows)} probe instances → {out_csv}")

def main():
    ap = argparse.ArgumentParser(description="Build a Packmol mixture around receptor for MixMD-style runs.")
    ap.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb")
    ap.add_argument("--box-mode", choices=["auto","fixed"], default="auto",
                    help="auto: size from receptor bounding box + padding; fixed: use --box-size-nm.")
    ap.add_argument("--box-size-nm", type=float, default=7.0,
                    help="Cubic box edge length if --box-mode fixed.")
    ap.add_argument("--padding-nm", type=float, default=1.2,
                    help="Padding added to receptor span if --box-mode auto.")
    ap.add_argument("--probe-total", type=int, default=200)
    ap.add_argument("--probe-fractions", type=str,
                    default="ipa=0.40,acn=0.20,imd=0.15,aceam=0.15,phol=0.10")
    ap.add_argument("--tolerance-A", type=float, default=2.0,
                    help="Packmol minimum distance in Å (applies between *any* atoms).")
    ap.add_argument("--output-prefix", type=str, default="mixmd")
    args = ap.parse_args()

    receptor_pdb = Path(args.input_receptor)
    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")

    build = Path("build"); tmp = build / "tmp_probes"
    build.mkdir(exist_ok=True); tmp.mkdir(parents=True, exist_ok=True)

    # Fractions → integer counts
    fr = parse_kv_fracs(args.probe_fractions)
    counts = {k: int(round(args.probe_total*fr[k])) for k in fr}
    # Ensure sum exactly equals probe_total
    diff = args.probe_total - sum(counts.values())
    if diff != 0:
        # adjust the largest fraction
        kmax = max(counts, key=lambda k: fr[k])
        counts[kmax] += diff

    # Determine box
    if args.box_mode == "auto":
        box = auto_box_from_receptor(receptor_pdb, args.padding_nm)
    else:
        # manual: side length in nm expected from --box-size-nm
        box = Box(side_nm=float(args.box_size_nm))
    print(f"[3c] Box: {box.L:.2f} nm (lo={box.lo:.2f}, hi={box.hi:.2f})")

    # Templates
    probe_pdbs = write_probe_templates(tmp)

    # Packmol input/output
    tag = args.output_prefix
    root = f"{receptor_pdb.stem}_{tag}"
    out_pdb = build / f"{root}.pdb"
    inp     = build / f"{root}.packmol.inp"
    # NOTE: pass a numeric side length (nm) here
    write_packmol_input(
        inp,
        receptor_pdb,
        probe_pdbs,
        counts,
        box.side_nm,          # <<< this is now guaranteed to exist
        args.tolerance_A,
        out_pdb
    )

    print(f"[3c] Running: packmol < {inp}")
    run_packmol(inp)

    # Centroids CSV for script 5
    out_csv = build / f"{receptor_pdb.stem}_{args.output-prefix}_placements.csv"
    probe_resnames = [PROBES[k][0] for k in fr.keys()]
    parse_probe_centroids(out_pdb, probe_resnames, out_csv)

    print("[3c] Done.")
    print(f"  • mixture pdb : {out_pdb}")
    print(f"  • placements  : {out_csv}")
    print(f"  • input file  : {inp}")

if __name__ == "__main__":
    main()
