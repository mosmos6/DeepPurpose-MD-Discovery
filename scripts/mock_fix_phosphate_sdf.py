#!/usr/bin/env python3
# One-off playground: ensure phosphate is represented as P=O to neutral O,
# and single bonds to O(-) in SDF. Then test OpenFF can read it.
#
# Usage:
#   conda run -n deeppurpose-md-env python scripts/mock_fix_phosphate_sdf.py \
#       --in fixed_UMP.pdb --out UMP_patched.sdf
#
# This does NOT modify your main step 4/5. It's only to discover a working recipe.

from __future__ import annotations
from pathlib import Path
import argparse, subprocess, sys

# ---- helpers for SDF parsing/writing (text-level; no RDKit required to fix) ----

def _read_sdf_blocks(sdf_path: Path):
    txt = sdf_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(txt) < 4:
        raise RuntimeError("SDF too short / missing header lines.")
    name, prog, comment = txt[0], txt[1], txt[2]
    counts = txt[3]
    # MDL counts line uses 3-char right-aligned integers
    try:
        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])
    except Exception:
        # robust fallback
        parts = counts.split()
        n_atoms, n_bonds = int(parts[0]), int(parts[1])

    atom_start = 4
    bond_start = atom_start + n_atoms
    atom_lines = txt[atom_start:bond_start]
    bond_lines = txt[bond_start:bond_start + n_bonds]
    tail_lines = txt[bond_start + n_bonds:]
    return (name, prog, comment, counts, atom_lines, bond_lines, tail_lines)

def _write_sdf_blocks(out_path: Path, name, prog, comment, counts, atom_lines, bond_lines, tail_lines):
    out = [name, prog, comment, counts] + atom_lines + bond_lines + tail_lines
    out_path.write_text("\n".join(out) + ("\n" if not out[-1].endswith("\n") else ""), encoding="utf-8")

def _atom_symbol(atom_line: str) -> str:
    # Columns 31-34 (0-based slice [31:34]) hold the atom symbol in V2000
    return atom_line[31:34].strip()

def _atom_charge_from_block(atom_line: str) -> int:
    # Atom block "charge code" is in columns 36-39 (0-based [36:39]); map code→formal charge
    # MDL mapping: 0→0, 1→+3, 2→+2, 3→+1, 4→0, 5→-1, 6→-2, 7→-3
    code_map = {0:0, 1:+3, 2:+2, 3:+1, 4:0, 5:-1, 6:-2, 7:-3}
    s = atom_line[36:39].strip()
    if not s:
        return 0
    try:
        return code_map.get(int(s), 0)
    except Exception:
        return 0

def _parse_m_chg(tail_lines):
    """Return dict {atom_idx: formal_charge} from any M  CHG lines (1-based indices)."""
    charges = {}
    for ln in tail_lines:
        if ln.startswith("M  CHG"):
            parts = ln.split()
            try:
                n = int(parts[2])
                for k in range(n):
                    ai = int(parts[3 + 2*k])      # 1-based
                    q  = int(parts[4 + 2*k])
                    charges[ai] = q
            except Exception:
                pass
    return charges

def _parse_bond_line(bline: str):
    # 3-digit right-aligned integer fields
    a = int(bline[0:3]); b = int(bline[3:6]); typ = int(bline[6:9]); stereo = int(bline[9:12])
    return a, b, typ, stereo

def _format_bond_line(old: str, a: int, b: int, typ: int, stereo: int):
    # preserve any trailing columns from old (beyond stereo)
    rest = old[12:]
    return f"{a:>3}{b:>3}{typ:>3}{stereo:>3}{rest}"

def _build_adjacency(n_atoms: int, bond_lines):
    adj = {i: [] for i in range(1, n_atoms+1)}
    for i, bl in enumerate(bond_lines):
        a, b, typ, stereo = _parse_bond_line(bl)
        adj[a].append((b, i))
        adj[b].append((a, i))
    return adj

def _effective_charges(atom_lines, tail_lines):
    """Combine atom-block charge code and M  CHG, with M  CHG taking precedence."""
    c_block = {i+1: _atom_charge_from_block(al) for i, al in enumerate(atom_lines)}
    c_mchg = _parse_m_chg(tail_lines)
    out = {}
    for i in range(1, len(atom_lines)+1):
        out[i] = c_mchg.get(i, c_block.get(i, 0))
    return out

def _fix_phosphates(atom_lines, bond_lines, tail_lines):
    """Edit bond orders so: each P has exactly one P=O to a neutral O, and all O(-) are single to P."""
    n_atoms = len(atom_lines)
    adj = _build_adjacency(n_atoms, bond_lines)
    charges = _effective_charges(atom_lines, tail_lines)

    # Gather symbols
    symbols = {i+1: _atom_symbol(al) for i, al in enumerate(atom_lines)}

    # Walk all P atoms
    for p_idx in range(1, n_atoms+1):
        if symbols[p_idx] != "P":
            continue
        # O neighbors of this P
        o_neighbors = [(nbr, bl_i) for (nbr, bl_i) in adj[p_idx] if symbols[nbr] == "O"]
        if not o_neighbors:
            continue

        # Identify which O are negative
        o_minus = {nbr for (nbr, _) in o_neighbors if charges.get(nbr, 0) < 0}
        o_neutral = [nbr for (nbr, _) in o_neighbors if nbr not in o_minus]

        # Choose which O should be P=O:
        # Prefer the O that is already double-bonded AND neutral; otherwise pick the first neutral.
        current_double = None
        for (nbr, bl_i) in o_neighbors:
            a, b, typ, stereo = _parse_bond_line(bond_lines[bl_i])
            if typ == 2:  # double
                current_double = nbr
                break
        if current_double is not None and current_double in o_neutral:
            double_target = current_double
        else:
            double_target = o_neutral[0] if o_neutral else None

        # Apply edits: all O(-) must be single; exactly one neutral O (double_target) is double; others single
        for (nbr, bl_i) in o_neighbors:
            a, b, typ, stereo = _parse_bond_line(bond_lines[bl_i])
            # normalize stereo to 0 for these bonds
            stereo = 0
            if nbr in o_minus:
                new_typ = 1
            elif double_target is not None and nbr == double_target:
                new_typ = 2
            else:
                new_typ = 1
            bond_lines[bl_i] = _format_bond_line(bond_lines[bl_i], a, b, new_typ, stereo)

    return bond_lines

# ---- main workflow ----

def main():
    ap = argparse.ArgumentParser(description="Mockup: make an SDF from PDB, patch phosphate valences, then test OpenFF reading.")
    ap.add_argument("--in", dest="inpdb", required=True, help="Input fixed ligand PDB (e.g., fixed_UMP.pdb)")
    ap.add_argument("--out", dest="outsdf", required=True, help="Output SDF path to write (patched)")
    ap.add_argument("--keep-initial", action="store_true", help="Also keep the initial OBabel SDF next to the patched one")
    args = ap.parse_args()

    inpdb = Path(args.inpdb).resolve()
    outsdf = Path(args.outsdf).resolve()
    work_sdf = outsdf.with_suffix(".initial.sdf")

    # 1) Produce initial SDF exactly like your single‑ligand pipeline
    cmd = ["obabel", inpdb, "-O", work_sdf]  # Path objects are fine; str() not required
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("❌ 'obabel' not found in PATH.", file=sys.stderr)
        sys.exit(2)
    except subprocess.CalledProcessError as e:
        print("❌ Open Babel failed to convert PDB→SDF.", file=sys.stderr)
        print((e.stderr or b"").decode(errors="ignore"), file=sys.stderr)
        sys.exit(2)

    # 2) Patch bond orders around phosphate(s) to be consistent with O(-) charges
    name, prog, comment, counts, atom_lines, bond_lines, tail_lines = _read_sdf_blocks(work_sdf)
    bond_lines = _fix_phosphates(atom_lines, bond_lines, tail_lines)
    _write_sdf_blocks(outsdf, name, prog, comment, counts, atom_lines, bond_lines, tail_lines)

    # 3) Test OpenFF can read the patched SDF
    try:
        from openff.toolkit.topology import Molecule
        m = Molecule.from_file(str(outsdf))
        if m.n_conformers == 0:
            print("ℹ️ OpenFF read the molecule but it had no conformers; adding one (for sanity).")
            m.generate_conformers(n_conformers=1)
        print(f"✅ OpenFF successfully read '{outsdf.name}'. Atoms: {m.n_atoms}, Bonds: {m.n_bonds}.")
    except Exception as e:
        print(f"❌ OpenFF could not read the patched SDF ({outsdf.name}).", file=sys.stderr)
        print(f"   Error: {e}", file=sys.stderr)
        sys.exit(1)

    # 4) Clean up initial SDF unless requested
    if not args.keep_initial and work_sdf.exists():
        try:
            work_sdf.unlink()
        except Exception:
            pass

if __name__ == "__main__":
    main()
