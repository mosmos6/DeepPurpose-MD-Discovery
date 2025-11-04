#!/usr/bin/env python3
# Mockup fixer: PDB -> SDF via Open Babel, then repair phosphate O(-) charges (V2000),
# then verify with OpenFF. No RDKit imports; textual SDF patch only.
#
# Usage examples:
#   python mock_fix_phosphate_sdf.py --in fixed_UMP.pdb --out UMP.sdf
#   python mock_fix_phosphate_sdf.py --in UMP.sdf --out UMP_fixed.sdf --skip-obabel
#
# What it does:
#   - If input is PDB and --skip-obabel is NOT set: obabel <pdb> -> temp.sdf (V2000)
#   - Parse SDF bonds to find P centers; for each P, find O neighbors
#       * O with degree==1 (only bonded to that P) are set to -1
#       * O with bond order 2 to P (P=O) is set to 0
#       * bridging O (degree>1) set to 0
#   - Writes a clean V2000 M  CHG line and sets atom-block charge code (5 for -1, 0 for neutral)
#   - Verifies OpenFF can read the final SDF
#
# Requirements in the environment:
#   - Open Babel CLI:   obabel
#   - OpenFF Toolkit:   openff-toolkit

from __future__ import annotations
from pathlib import Path
import argparse
import subprocess
import sys
import shutil

# OpenFF just for validation (no RDKit import here)
from openff.toolkit.topology import Molecule
from openff.toolkit.utils.exceptions import MoleculeParseError


def run_obabel_pdb_to_sdf(pdb: Path, sdf_out: Path) -> None:
    """Convert PDB -> SDF (V2000) using Open Babel, mirroring the single‑ligand pipeline."""
    if shutil.which("obabel") is None:
        sys.exit("ERROR: Open Babel ('obabel') not found in PATH. Install it in your env.")
    # Simple, conservative conversion (no V3000 by default; matches your older flow)
    cmd = ["obabel", str(pdb), "-O", str(sdf_out)]
    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        msg = (e.stderr or e.stdout or str(e)).strip()
        sys.exit(f"ERROR: obabel failed on {pdb.name} -> {sdf_out.name}\n{msg}")


def _parse_v2000_counts_line(line: str) -> tuple[int, int]:
    """Return (natoms, nbonds) from the V2000 counts line."""
    # V2000 counts line: columns are fixed-width; but splitting works for OB output.
    parts = line.split()
    if len(parts) < 2:
        raise ValueError("Malformed V2000 counts line")
    try:
        na = int(parts[0]); nb = int(parts[1])
    except Exception as e:
        raise ValueError(f"Malformed counts values: {line!r}") from e
    return na, nb


def _read_sdf_v2000(path: Path):
    """
    Read single‑molecule SDF (V2000) into:
      atoms: list of dict(idx=1-based, symbol, chg_code[int])
      bonds: list of dict(a1,a2,order)
      header_lines: first 3 lines (title, program, comment)
      counts_line: the counts line (string)
      atom_lines, bond_lines: original text blocks (lists of strings)
      tail_lines: lines after bond block (incl. M  CHG if present, M  END, $$$$)
    """
    txt = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=False)

    if len(txt) < 4:
        raise ValueError("SDF too short / not V2000")

    header_lines = txt[0:3]
    counts_line = txt[3]
    natoms, nbonds = _parse_v2000_counts_line(counts_line)

    atom_start = 4
    atom_end = atom_start + natoms
    bond_start = atom_end
    bond_end = bond_start + nbonds

    atom_lines = txt[atom_start:atom_end]
    bond_lines = txt[bond_start:bond_end]
    tail_lines = txt[bond_end:]

    # Parse atoms
    atoms = []
    for i, ln in enumerate(atom_lines, start=1):
        # x y z symbol massDiff charge ...  -> symbol at col 31-34 in strict format; OB outputs space-delimited
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed atom line {i}: {ln!r}")
        symbol = parts[3]
        chg_code = 0
        if len(parts) >= 6:
            # parts[5] is the "charge code" for V2000 if OB wrote it
            # (codes: 0 none, 1:+3, 2:+2, 3:+1, 4:doublet, 5:-1, 6:-2, 7:-3)
            try:
                chg_code = int(parts[5])
            except Exception:
                chg_code = 0
        atoms.append({"idx": i, "symbol": symbol, "chg_code": chg_code, "raw": ln})

    # Parse bonds
    bonds = []
    for ln in bond_lines:
        parts = ln.split()
        if len(parts) < 3:
            raise ValueError(f"Malformed bond line: {ln!r}")
        a1 = int(parts[0]); a2 = int(parts[1]); order = int(parts[2])
        bonds.append({"a1": a1, "a2": a2, "order": order, "raw": ln})

    return {
        "atoms": atoms,
        "bonds": bonds,
        "header_lines": header_lines,
        "counts_line": counts_line,
        "atom_lines": atom_lines,
        "bond_lines": bond_lines,
        "tail_lines": tail_lines,
    }


def _write_sdf_v2000(path: Path, data, new_MCHG: dict[int, int]) -> None:
    """
    Write V2000 SDF with:
      - atom-block 'charge code' updated to match new_MCHG (-1 -> code 5, 0 -> 0)
      - a single normalized 'M  CHG' line reflecting new_MCHG (only nonzero charges)
    """
    atoms = data["atoms"]
    bonds = data["bonds"]
    header_lines = data["header_lines"]
    counts_line = data["counts_line"]
    tail_lines = data["tail_lines"]

    # Rebuild atom lines with corrected charge code
    new_atom_lines = []
    for atom in atoms:
        ln = atom["raw"]
        code = 5 if new_MCHG.get(atom["idx"], 0) == -1 else 0
        # Rebuild the 6th numeric field (charge code) if present; else append
        parts = ln.split()
        if len(parts) >= 6:
            parts[5] = str(code)
            # Recreate an OB-like spacing; accept looser whitespace (OpenFF is fine)
            new_ln = " ".join(parts)
        else:
            new_ln = ln  # Extremely unlikely with OB output
        new_atom_lines.append(new_ln)

    # Keep original bond lines
    new_bond_lines = [b["raw"] for b in bonds]

    # Strip any existing M  CHG lines in tail; keep others (M  END, $$$$, etc.)
    kept_tail = [ln for ln in tail_lines if not ln.startswith("M  CHG")]

    # Build new M  CHG (only nonzero)
    entries = [(idx, q) for idx, q in sorted(new_MCHG.items()) if q != 0]
    if entries:
        # M  CHG  n  a1 q1  a2 q2 ...
        parts = ["M  CHG", f"{len(entries):>3}"]
        for (aidx, q) in entries:
            parts.append(f"{aidx:>4}")
            parts.append(f"{q:>4}")
        mchg_line = " ".join(parts)
        # Insert before M  END if present; else append
        inserted = False
        new_tail = []
        for ln in kept_tail:
            if (not inserted) and ln.startswith("M  END"):
                new_tail.append(mchg_line)
                inserted = True
            new_tail.append(ln)
        if not inserted:
            new_tail.append(mchg_line)
        kept_tail = new_tail

    out_lines = []
    out_lines.extend(header_lines)
    out_lines.append(counts_line)
    out_lines.extend(new_atom_lines)
    out_lines.extend(new_bond_lines)
    out_lines.extend(kept_tail)
    if not any(ln.strip() == "$$$$" for ln in kept_tail):
        out_lines.append("$$$$")

    path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def _build_neighbors(n_atoms: int, bonds: list[dict]) -> list[list[tuple[int,int]]]:
    """Return neighbor list: for each atom idx (1-based), list of (nbr_idx, order)."""
    nbrs = [[] for _ in range(n_atoms + 1)]
    for b in bonds:
        a1, a2, order = b["a1"], b["a2"], b["order"]
        nbrs[a1].append((a2, order))
        nbrs[a2].append((a1, order))
    return nbrs


def _repair_phosphate_charges_v2000(data) -> dict[int, int]:
    """
    Decide charges purely from topology:
      For each P:
        - O neighbors with degree==1 (only bonded to P): set -1
        - O neighbor with P=O (order==2): set 0
        - O neighbors with degree>1 (bridging): set 0
      Everything else: 0
    Return a dict {atom_idx: formal_charge}.
    """
    atoms = data["atoms"]
    bonds = data["bonds"]
    na = len(atoms)
    nbrs = _build_neighbors(na, bonds)

    # Start neutral
    chg = {i: 0 for i in range(1, na + 1)}

    # Find P centers
    p_indices = [a["idx"] for a in atoms if a["symbol"].upper() == "P"]
    if not p_indices:
        # Nothing to do
        return chg

    for p in p_indices:
        # Look at P's O neighbors
        onbr = [(j, order) for (j, order) in nbrs[p] if atoms[j-1]["symbol"].upper() == "O"]

        # Classify O’s
        singly_bound = []
        double_bonded = []
        bridging = []

        for (o_idx, order) in onbr:
            deg = len(nbrs[o_idx])
            if order == 2:
                double_bonded.append(o_idx)
            elif deg == 1:
                singly_bound.append(o_idx)
            else:
                bridging.append(o_idx)

        # Heuristic for phosphate: 1 double-bonded O (neutral), 2 singly-bound O (-1),
        # 1 bridging O to sugar (neutral). If counts mismatch, we still apply the rule we can.
        for o in double_bonded:
            chg[o] = 0
        for o in bridging:
            chg[o] = 0
        for o in singly_bound:
            chg[o] = -1  # <- THE FIX

        # Optional sanity: if we didn’t find 2 singly-bound O, warn (but still proceed)
        if len(singly_bound) != 2:
            sys.stderr.write(
                f"⚠️  Warning: P atom {p} has {len(singly_bound)} singly-bound O "
                f"(expected 2). Will still write charges based on this classification.\n"
            )

    return chg


def _openff_can_read(sdf: Path) -> tuple[bool, str]:
    try:
        _ = Molecule.from_file(str(sdf))
        return True, "ok"
    except Exception as e:
        return False, str(e)


def main():
    ap = argparse.ArgumentParser(description="Mockup: make SDF for UMP and fix phosphate O(-) placement")
    ap.add_argument("--in", dest="inp", required=True, help="Input file: PDB (typical) or SDF")
    ap.add_argument("--out", dest="out", required=True, help="Output SDF path to write")
    ap.add_argument("--skip-obabel", action="store_true", help="Assume input is already SDF; skip PDB->SDF conversion")
    ap.add_argument("--keep-temp", action="store_true", help="Keep intermediate SDF before patching")
    args = ap.parse_args()

    in_path = Path(args.inp).resolve()
    out_sdf = Path(args.out).resolve()

    if not in_path.exists():
        sys.exit(f"ERROR: Input file not found: {in_path}")

    # 1) Produce an SDF (V2000) via OBabel unless told to skip
    if args.skip_obabel:
        if in_path.suffix.lower() != ".sdf":
            sys.exit("ERROR: --skip-obabel was set but input is not an SDF.")
        tmp_sdf = in_path
    else:
        if in_path.suffix.lower() != ".pdb":
            sys.exit("ERROR: Input must be a PDB unless --skip-obabel is passed.")
        tmp_sdf = out_sdf.with_suffix(".prepatch.sdf")
        run_obabel_pdb_to_sdf(in_path, tmp_sdf)

    # 2) Read SDF (V2000), repair phosphate charges
    data = _read_sdf_v2000(tmp_sdf)
    new_chg = _repair_phosphate_charges_v2000(data)
    _write_sdf_v2000(out_sdf, data, new_chg)

    # 3) Validate with OpenFF
    ok, msg = _openff_can_read(out_sdf)
    if not ok:
        sys.stderr.write(f"❌ OpenFF could not read the patched SDF ({out_sdf.name}).\n")
        sys.stderr.write(f"   Error: {msg}\n")
        sys.exit(1)

    # Cleanup temp
    if (not args.keep_temp) and tmp_sdf != out_sdf and tmp_sdf.exists():
        try: tmp_sdf.unlink()
        except Exception: pass

    # Report what we did (which O’s got -1)
    neg_idx = [i for i, q in new_chg.items() if q == -1]
    print(f"✅ Wrote {out_sdf.name} and verified with OpenFF.")
    if neg_idx:
        print(f"   Set formal charge −1 on O atoms (by index): {', '.join(map(str, neg_idx))}")
    else:
        print("   No phosphate O(−) atoms detected; wrote neutral charges.")

if __name__ == "__main__":
    main()
