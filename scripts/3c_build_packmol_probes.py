#!/usr/bin/env python3
# scripts/3c_build_packmol_probes.py
# Universal Packmol driver for MixMD-style placement (protein/RNA targets, incl. 40S)
# No OpenMM/OpenFF/JAX here — RDKit + stdlib only.
# Updated: 2025-10-15

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass

# ---- RDKit (templates only) --------------------------------------------------
from rdkit import Chem
from rdkit.Chem import AllChem

# ---- Probes: key -> (RESNAME, SMILES) ---------------------------------------
PROBES: Dict[str, Tuple[str, str]] = {
    "ipa":   ("IPA",  "CC(C)O"),        # isopropanol
    "acn":   ("ACN",  "CC#N"),          # acetonitrile
    "imd":   ("IMD",  "c1[nH]cnc1"),    # imidazole (neutral, N1 protonated)
    "aceam": ("ACEA", "CC(=O)N"),       # acetamide
    "phol":  ("PHOL", "c1ccc(cc1)O"),   # phenol
    "acoh":  ("ACOH", "CC(=O)O"),       # acetic acid
}

# ---- Single box type (nm) ----------------------------------------------------
@dataclass(frozen=True)
class CubicBox:
    """Simple cubic box in nm."""
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

# ---- Helpers ----------------------------------------------------------------
def parse_kv_fracs(s: str) -> Dict[str, float]:
    """Parse 'ipa=0.40,acn=0.20,...' -> normalized fractions."""
    fr: Dict[str, float] = {}
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
    """
    Build a tiny PDB block for one small molecule, centered at its centroid,
    with residue name = resname. Coordinates: Å.
    """
    mol = Chem.AddHs(mol)
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    coords = []
    cx = cy = cz = 0.0
    for i in range(n):
        pos = conf.GetAtomPosition(i)  # RDKit Point3D (Å)
        x, y, z = float(pos.x), float(pos.y), float(pos.z)
        coords.append((x, y, z))
        cx += x; cy += y; cz += z
    cx /= n; cy /= n; cz /= n

    lines = []
    serial = 1
    for i, (x, y, z) in enumerate(coords):
        atom = mol.GetAtomWithIdx(i)
        elem = atom.GetSymbol().rjust(2)
        name = (elem + "  ").ljust(4)[:4]
        x0 = x - cx; y0 = y - cy; z0 = z - cz
        # HETATM serial, atom name, resname, chain/resid, x y z, occ, b, element
        lines.append(f"HETATM{serial:5d} {name} {resname:>3s} Z{1:4d}    "
                     f"{x0:8.3f}{y0:8.3f}{z0:8.3f}  1.00  0.00          {elem:>2s}")
        serial += 1
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"

def write_probe_templates(out_dir: Path) -> Dict[str, Path]:
    """Create one centered PDB per probe in out_dir and return key->path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Path] = {}
    for key, (resname, smi) in PROBES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise RuntimeError(f"SMILES parse failed for {key}: {smi}")
        pdb_block = rdkit_centered_pdb(mol, resname)
        p = out_dir / f"{resname}.pdb"
        p.write_text(pdb_block, encoding="utf-8")
        out[key] = p
    return out

def auto_box_from_receptor(receptor_pdb: Path, padding_nm: float) -> CubicBox:
    """
    Compute a cubic box that encloses the receptor + padding_nm (nm).
    PDB coords are Å; convert to nm.
    """
    coords_A: List[Tuple[float, float, float]] = []
    with open(receptor_pdb, "r") as fh:
        for line in fh:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 54:
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    coords_A.append((x, y, z))
                except ValueError:
                    pass
    if not coords_A:
        raise ValueError("No ATOM/HETATM coordinates found in receptor PDB.")

    import numpy as _np
    xyz_nm = _np.asarray(coords_A, dtype=_np.float64) * 0.1  # Å → nm
    mins = xyz_nm.min(axis=0)
    maxs = xyz_nm.max(axis=0)
    span_nm = float((maxs - mins).max())
    side_nm = span_nm + 2.0 * float(padding_nm)

    half = 0.5 * side_nm
    print(f"[3c] Box: {side_nm:.2f} nm (lo={-half:.2f}, hi={half:.2f})")
    return CubicBox(side_nm=side_nm)

def write_packmol_input(inp_path: Path,
                        receptor_pdb: Path,
                        probe_pdbs: Dict[str, Path],
                        counts: Dict[str, int],
                        box_side_nm: float,
                        tol_A: float,
                        out_pdb: Path) -> None:
    """
    Emit a Packmol input file with:
      - receptor fixed at the origin,
      - probes uniformly inside the full cubic box.
    All distances passed to Packmol are in Å.
    """
    side_A = float(box_side_nm) * 10.0
    half_A = 0.5 * side_A
    lo, hi = -half_A, half_A

    lines: List[str] = []
    lines.append(f"tolerance {float(tol_A):.3f}")
    lines.append("filetype pdb")
    lines.append("add_box_sides")
    lines.append(f"output {out_pdb.as_posix()}")

    # Receptor fixed at center
    lines.append(f"structure {receptor_pdb.as_posix()}")
    lines.append("  number 1")
    lines.append("  center")
    lines.append("  fixed 0. 0. 0. 0. 0. 0.")
    lines.append("end structure")

    # Probes across the full box
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

def run_packmol(input_file: Path) -> None:
    """
    Run Packmol as: packmol < input_file, with a *seekable* stdin handle
    to avoid Fortran REWIND issues.
    """
    exe = shutil.which("packmol")
    if exe is None:
        raise RuntimeError("Packmol executable not found in PATH.")
    print(f"[3c] Running: {exe} < {input_file}")

    with open(input_file, "rb") as fin:
        proc = subprocess.run(
            [exe],
            stdin=fin,  # seekable handle
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )

    if proc.returncode != 0:
        sys.stderr.write(proc.stdout.decode("utf-8", errors="ignore"))
        sys.stderr.write(proc.stderr.decode("utf-8", errors="ignore"))
        raise RuntimeError("[3c] Packmol failed. See console above.")
    print("[3c] Packmol completed.")

def parse_probe_centroids(packmol_pdb: Path, probe_resnames: List[str], out_csv: Path) -> None:
    """Compute centroids (nm) per (resname,resid) from the Packmol mixture PDB."""
    groups: Dict[Tuple[str, int], List[Tuple[float, float, float]]] = {}
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

    rows: List[Tuple[str, int, float, float, float]] = []
    for (resname, resid), xyzs in groups.items():
        n = len(xyzs)
        cx = sum(p[0] for p in xyzs)/n * 0.1  # nm
        cy = sum(p[1] for p in xyzs)/n * 0.1
        cz = sum(p[2] for p in xyzs)/n * 0.1
        rows.append((resname, resid, cx, cy, cz))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["resname", "resid", "x_nm", "y_nm", "z_nm"])
        w.writerows(rows)

    print(f"[3c] Wrote centroids for {len(rows)} probe instances → {out_csv.as_posix()}")

# ---- CLI --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Build a Packmol mixture around receptor for MixMD-style runs.")
    ap.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb")
    ap.add_argument("--box-mode", choices=["auto", "fixed"], default="auto",
                    help="auto = fit receptor bbox + padding; fixed = use --box-size-nm.")
    ap.add_argument("--box-size-nm", type=float, default=7.0,
                    help="Cubic box edge length if --box-mode fixed.")
    ap.add_argument("--padding-nm", type=float, default=1.2,
                    help="Extra margin added to receptor span if --box-mode auto.")
    ap.add_argument("--probe-total", type=int, default=200)
    ap.add_argument("--probe-fractions", type=str,
                    default="ipa=0.40,acn=0.20,imd=0.15,aceam=0.15,phol=0.10")
    ap.add_argument("--tolerance-A", type=float, default=2.0,
                    help="Packmol minimum atom–atom distance in Å (2.0–2.2 typical).")
    ap.add_argument("--output-prefix", type=str, default="mixmd")
    args = ap.parse_args()

    receptor_pdb = Path(args.input_receptor)
    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")

    build = Path("build")
    tmp   = build / "tmp_probes"
    build.mkdir(exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)

    # Fractions → integer counts (exactly sums to probe_total)
    fr = parse_kv_fracs(args.probe_fractions)
    counts = {k: int(round(args.probe_total * fr[k])) for k in fr}
    diff = args.probe_total - sum(counts.values())
    if diff != 0:
        kmax = max(counts, key=lambda k: fr[k])
        counts[kmax] += diff

    # Box
    if args.box_mode == "auto":
        box = auto_box_from_receptor(receptor_pdb, args.padding_nm)
    else:
        box = CubicBox(side_nm=float(args.box_size_nm))
        half = box.half_nm
        print(f"[3c] Box: {box.side_nm:.2f} nm (lo={-half:.2f}, hi={half:.2f})")

    # Probe templates
    probe_pdbs = write_probe_templates(tmp)

    # Packmol I/O
    tag  = args.output_prefix
    root = f"{receptor_pdb.stem}_{tag}"
    out_pdb = build / f"{root}.pdb"
    inp     = build / f"{root}.packmol.inp"

    write_packmol_input(
        inp,
        receptor_pdb,
        probe_pdbs,
        counts,
        box.side_nm,            # numeric side length in nm
        args.tolerance_A,
        out_pdb
    )

    print(f"[3c] Running: packmol < {inp}")
    run_packmol(inp)

    # Centroids CSV for downstream scripts
    out_csv = build / f"{receptor_pdb.stem}_{tag}_placements.csv"
    probe_resnames = [PROBES[k][0] for k in fr.keys()]
    parse_probe_centroids(out_pdb, probe_resnames, out_csv)

    print("[3c] Done.")
    print(f"  • mixture pdb : {out_pdb}")
    print(f"  • placements  : {out_csv}")
    print(f"  • input file  : {inp}")

if __name__ == "__main__":
    main()
