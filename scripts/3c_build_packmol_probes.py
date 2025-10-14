#!/usr/bin/env python3
# scripts/3c_build_packmol_probes.py
# Universal Packmol mixture builder (target-agnostic)
# Creates a cosolvent PDB around the receptor so coordinates align with your MD.

import argparse, os, shutil, subprocess, sys
from pathlib import Path

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openmm.app import PDBFile, Modeller
from openmm import Vec3
from openmm.unit import angstrom, nanometer

PROBES = {
    # key      resname  smiles
    "ipa":   ("IPA",  "CC(C)O"),            # isopropanol
    "acn":   ("ACN",  "CC#N"),              # acetonitrile
    "imd":   ("IMD",  "c1[nH]cnc1"),        # imidazole (N1 protonated)
    "aceam": ("ACEA", "CC(=O)N"),           # acetamide
    "phol":  ("PHOL", "c1ccc(cc1)O"),       # phenol
    "acoh":  ("ACOH", "CC(=O)O"),           # acetic acid (optional)
}

def parse_fractions(s: str):
    # accepts "ipa=0.40,acn=0.20,imd=0.15,aceam=0.15,phol=0.10"
    parts = [p.strip() for p in s.split(",") if p.strip()]
    kv = {}
    for p in parts:
        if "=" not in p:
            raise ValueError(f"Bad token: {p!r} (use key=val)")
        k, v = p.split("=", 1)
        k = k.strip().lower()
        if k not in PROBES:
            raise ValueError(f"Unknown key {k!r}. Allowed: {list(PROBES.keys())}")
        kv[k] = float(v.strip())
    total = sum(kv.values())
    if total <= 0:
        raise ValueError("Probe fractions must sum to > 0")
    # normalize so minor rounding doesn’t matter
    for k in kv:
        kv[k] /= total
    return kv

def write_probe_pdbs(out_dir: Path):
    """Write one PDB template per probe with proper residue name."""
    out_dir.mkdir(parents=True, exist_ok=True)
    probe_pdb_paths = {}
    for key, (resname, smiles) in PROBES.items():
        off = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = resname
        top = OFFTopology.from_molecules([off]).to_openmm()
        # Center on own centroid (Å -> nm inside write)
        conf = off.conformers[0]  # quantity (Å)
        # Convert Å to nm and center
        coords_nm = []
        cx = sum([v.x for v in conf.value_in_unit(angstrom)]) / len(conf)
        cy = sum([v.y for v in conf.value_in_unit(angstrom)]) / len(conf)
        cz = sum([v.z for v in conf.value_in_unit(angstrom)]) / len(conf)
        for v in conf.value_in_unit(angstrom):
            coords_nm.append(Vec3((v.x - cx)*0.1, (v.y - cy)*0.1, (v.z - cz)*0.1))  # nm
        mod = Modeller(top, coords_nm)
        # rename the single residue to resname
        for res in mod.topology.residues():
            res.name = resname
        out_pdb = out_dir / f"{resname}.pdb"
        with open(out_pdb, "w") as f:
            PDBFile.writeFile(mod.topology, mod.positions, f, keepIds=True)
        probe_pdb_paths[key] = out_pdb
    return probe_pdb_paths

def build_packmol_input(
    receptor_pdb: Path,
    probe_pdbs: dict,
    counts: dict,
    box_size_nm: float,
    tolerance_ang: float,
    out_pdb: Path,
    inp_path: Path,
):
    """Compose a Packmol input that fixes the receptor and packs probes inside the same 0..L nm cube."""
    L_ang = box_size_nm * 10.0  # nm->Å; Packmol uses Å
    lines = []
    lines.append(f"tolerance {tolerance_ang:.3f}")
    lines.append("filetype pdb")
    lines.append(f"output {out_pdb.as_posix()}")
    lines.append("add_box_sides 0.0")  # keep exact box
    lines.append("")  # receptor fixed to preserve frame
    lines.append(f"struct {receptor_pdb.as_posix()}")
    lines.append("  number 1")
    lines.append("  fixed 0. 0. 0. 0. 0. 0.")
    lines.append("end struct")
    lines.append("")
    # probes
    for key, n in counts.items():
        if n <= 0:
            continue
        resname = PROBES[key][0]
        lines.append(f"struct {probe_pdbs[key].as_posix()}")
        lines.append(f"  number {int(n)}")
        lines.append(f"  inside box 0.0 0.0 0.0 {L_ang:.3f} {L_ang:.3f} {L_ang:.3f}")
        lines.append("end struct")
        lines.append("")
    inp_path.write_text("\n".join(lines))

def main():
    ap = argparse.ArgumentParser(description="Build a Packmol cosolvent mixture PDB (universal, target‑agnostic).")
    ap.add_argument("--input-receptor", default="receptor_cleaned.pdb")
    ap.add_argument("--target-name", default=None, help="Name used in outputs; default = stem of receptor filename")
    ap.add_argument("--box-size-nm", type=float, default=7.0, help="Cubic box edge used in MD and packing (nm)")
    ap.add_argument("--tolerance-ang", type=float, default=2.0, help="Packmol inter‑molecule tolerance (Å)")
    ap.add_argument("--probe-total", type=int, default=200, help="Total probes to place")
    ap.add_argument("--probe-fractions", default="ipa=0.40,acn=0.20,imd=0.15,aceam=0.15,phol=0.10",
                    help="Comma‑sep key=frac. Keys: ipa,acn,imd,aceam,phol,acoh")
    ap.add_argument("--packmol-bin", default="packmol", help="Packmol executable on PATH")
    ap.add_argument("--out-dir", default="build", help="Output directory")
    args = ap.parse_args()

    rec_path = Path(args.input_receptor)
    if not rec_path.exists():
        sys.exit(f"[3c] Receptor not found: {rec_path}")

    target = args.target_name or rec_path.stem  # universal — no hard‑coded target strings
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "packmol_inputs"
    tmp_dir.mkdir(exist_ok=True)

    # 1) make per‑probe PDB templates
    probe_pdbs = write_probe_pdbs(tmp_dir)

    # 2) decide counts from fractions
    fracs = parse_fractions(args.probe_fractions)
    counts = {k: int(round(args.probe_total * fracs.get(k, 0.0))) for k in PROBES.keys()}
    # last key absorbs rounding
    deficit = args.probe_total - sum(counts.values())
    if deficit != 0 and fracs:
        last = list(fracs.keys())[-1]
        counts[last] += deficit

    # 3) write Packmol input
    out_pdb = out_dir / f"{target}_mixmd.probes.pdb"
    inp_path = out_dir / f"{target}_mixmd.packmol.inp"
    build_packmol_input(rec_path, probe_pdbs, counts, args.box_size_nm, args.tolerance_ang, out_pdb, inp_path)

    # 4) run Packmol
    cmd = [args.packmol_bin, "<", inp_path.as_posix()]
    print(f"[3c] Running: {' '.join(cmd)}")
    # Use shell so "<" redirection works cross‑platform where Packmol expects stdin
    ret = subprocess.run(" ".join(cmd), shell=True)
    if ret.returncode != 0 or not out_pdb.exists():
        sys.exit("[3c] Packmol failed (see console).")
    print(f"[3c] ✅ Wrote mixture: {out_pdb.as_posix()}")
    print(f"[3c] Next: run 5 with --packmol-probes-pdb {out_pdb.as_posix()}")

if __name__ == "__main__":
    main()
