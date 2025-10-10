# scripts/mixmd_utils_5c.py
# Date 10/10/2025
# Utilities for MixMD-lite: parsing probe spec, building OFF molecules, and water->probe replacement.
# No "from openmm.unit import *" here to avoid shadowing builtins.

from __future__ import annotations
import builtins as _bi
import numpy as _np
import random
from typing import Dict, List, Tuple

import openmm as mm
from openmm import Vec3
from openmm.unit import angstrom, nanometer, Quantity

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# Minimal probe library (small, chemically diverse)
PROBES = {
    "ipa":   {"resname": "IPA",   "smiles": "CC(C)O"},       # isopropanol
    "acn":   {"resname": "ACN",   "smiles": "CC#N"},         # acetonitrile
    "imd":   {"resname": "IMD",   "smiles": "c1[nH]cnc1"},   # imidazole (neutral)
    "aceam": {"resname": "ACEA",  "smiles": "CC(=O)N"},      # acetamide
    "phol":  {"resname": "PHOL",  "smiles": "c1ccc(cc1)O"},  # phenol
    "acoh":  {"resname": "ACOH",  "smiles": "CC(=O)O"}       # acetic acid
}

WATER_NAMES = {"HOH", "WAT", "T3P", "TIP3", "SOL", "OH2"}

def parse_probe_fractions(spec: str) -> Dict[str, float]:
    """
    Parse "ipa=0.40,acn=0.20,..." or "ipa:0.40,..." into a normalized dict.
    Uses builtins.sum to avoid OpenMM unit math collisions.
    """
    kv: Dict[str, float] = {}
    for token in spec.split(","):
        t = token.strip()
        if not t:
            continue
        if "=" in t:
            k, v = t.split("=", 1)
        elif ":" in t:
            k, v = t.split(":", 1)
        else:
            raise ValueError(f"Bad token in --probe-fractions: '{t}' (use key=val or key:val)")
        k = k.strip()
        v = float(v.strip())
        if v < 0:
            raise ValueError(f"Negative fraction for '{k}'")
        kv[k] = v
    total = float(_bi.sum(kv.values()))  # builtins.sum
    if total <= 0.0:
        raise ValueError("Total probe fractions must be > 0")
    for k in list(kv.keys()):
        kv[k] = kv[k] / total
    return kv

def build_probe_molecules(keys: List[str]) -> Dict[str, dict]:
    """
    For each probe key, build an OpenFF Molecule (1 conformer), prepare an OpenMM Topology,
    and pre-center its conformer in *nanometers* around origin for easy translation.
    Returns dict[key] = {off, top, pos_nm (np.ndarray[n,3]), resname, natoms}
    """
    cache: Dict[str, dict] = {}
    for key in keys:
        if key not in PROBES:
            raise ValueError(f"Unknown probe key '{key}'. Allowed: {sorted(PROBES.keys())}")
        info = PROBES[key]
        off = Molecule.from_smiles(info["smiles"], allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = info["resname"]

        xyz_ang = _np.asarray(off.conformers[0].value_in_unit(angstrom))  # (n,3), plain floats
        centroid = xyz_ang.mean(axis=0, keepdims=True)
        xyz_nm = (xyz_ang - centroid) * 0.1  # Å -> nm, centered

        top = OFFTopology.from_molecules([off]).to_openmm()
        cache[key] = {
            "off": off,
            "top": top,
            "pos_nm": xyz_nm,
            "resname": info["resname"],
            "natoms": xyz_nm.shape[0],
        }
    return cache

def register_smirnoff_templates(forcefield, molecules: List[Molecule]) -> None:
    """Register a single SMIRNOFFTemplateGenerator for all molecules (ligand+probes)."""
    gen = SMIRNOFFTemplateGenerator(molecules=molecules)
    forcefield.registerTemplateGenerator(gen.generator)

def _count_atoms_in_residue(res) -> int:
    return _bi.sum(1 for _ in res.atoms())

def replace_waters_with_probes(
    modeller,
    probe_cache: Dict[str, dict],
    probe_fracs: Dict[str, float],
    total_count: int,
    seed: int = 24680,
    jitter_nm: float = 0.05,
) -> None:
    """
    Select random water residues and replace each with a probe (types distributed by probe_fracs).
    Operates *after* solvation. Positions are given as Quantity(list[Vec3], nanometer).
    """
    rng = random.Random(int(seed))

    waters = [res for res in modeller.topology.residues()
              if res.name in WATER_NAMES and _count_atoms_in_residue(res) >= 3]
    n_waters = len(waters)
    if n_waters == 0:
        print("⚠️ No recognizable waters (HOH/WAT/TIP3/T3P/SOL/OH2). Skipping MixMD placement.")
        return

    # Compute integer counts per probe (last bucket absorbs rounding)
    keys = list(probe_fracs.keys())
    counts = {}
    remaining = total_count
    for i, k in enumerate(keys):
        if i < len(keys) - 1:
            n = int(round(total_count * float(probe_fracs[k])))
            counts[k] = max(0, n)
            remaining -= n
        else:
            counts[k] = max(0, remaining)

    n_requested = int(_bi.sum(counts.values()))
    if n_requested == 0:
        print("⚠️ --probe-total resolved to 0 after rounding. Nothing to place.")
        return
    if n_requested > n_waters:
        raise ValueError(f"Requested {n_requested} probes but only {n_waters} waters are present.")

    rng.shuffle(waters)
    # Build a flat list of probe keys (e.g., ["ipa","ipa","acn",...]) and pair with first n_requested waters
    probe_keys_stream: List[str] = []
    for k, n in counts.items():
        probe_keys_stream.extend([k] * n)
    rng.shuffle(probe_keys_stream)

    to_delete = []
    to_add: List[Tuple] = []  # (topology, positions Quantity[nm])

    # We’ll translate each centered probe conformer onto the oxygen position of the chosen water
    for res, k in zip(waters[:n_requested], probe_keys_stream):
        # oxygen coordinate
        oxy = next(a for a in res.atoms() if a.element.symbol == "O")
        wpos = modeller.positions[oxy.index]  # Vec3 with units
        t = _np.array([float(wpos.x), float(wpos.y), float(wpos.z)])

        jitter = _np.array([rng.uniform(-jitter_nm, jitter_nm),
                            rng.uniform(-jitter_nm, jitter_nm),
                            rng.uniform(-jitter_nm, jitter_nm)])

        # place probe at water oxygen + jitter
        placed_nm = probe_cache[k]["pos_nm"] + (t + jitter)  # (natoms,3)
        vecs = [Vec3(float(x), float(y), float(z)) for (x, y, z) in placed_nm]
        pos_q = Quantity(vecs, nanometer)

        to_add.append((probe_cache[k]["top"], pos_q))
        to_delete.extend(list(res.atoms()))

    # Single delete, then add all probes
    modeller.delete(to_delete)
    for top, pos_q in to_add:
        modeller.add(top, pos_q)

    print(f"✅ MixMD-lite: placed {n_requested} probes "
          f"({', '.join([f'{k}:{counts[k]}' for k in keys if counts[k]>0])}).")
