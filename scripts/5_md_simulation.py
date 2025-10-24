# 5_md_simulation.py
# Author: Iori Mochizuki + collaborator
# Updated: 2025-10-24
# Description: OpenMM MD (protein/RNA/multi‑ligand) with optional MixMD-from-Packmol import.
# - MixMD: reads build/<receptor>_mixmd_placements.csv, crops to MD box, parametrizes probes (SMIRNOFF),
#   adds them before solvation, logs probe centroids over time (receptor-centered) to CSV.
# - Multi‑ligand: auto-detect ligand_*.sdf (or use --input-ligand), parametrizes & adds all ligands.
# - Robust to OpenMM unit-math shadowing of Python builtins (uses plain loops / builtins explicitly).
# - Writes both PDB snapshots and final mmCIF (avoids 5-digit PDB serial overflow).

import argparse, csv, glob, builtins as _py
from pathlib import Path
from sys import stdout

import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
import openmm.unit as u  # avoid wildcard; use u.nanometer, etc.

from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---------------------------
# Utilities
# ---------------------------

def _count_atoms(top: Topology) -> int:
    # Avoid OpenMM's unit-math sum() by using a plain loop
    n = 0
    for _ in top.atoms():
        n += 1
    return n

def _list_residue_names(top: Topology):
    names = []
    for r in top.residues():
        names.append((r.chain.id if hasattr(r, "chain") else "", (r.name or "").upper().strip()))
    return names

def _glob_ligands(default: str | None):
    paths = []
    if default:
        p = Path(default)
        if p.exists():
            paths.append(p)
    # auto-discover ligand_*.sdf
    for p in sorted(Path(".").glob("ligand_*.sdf")):
        if p not in paths:
            paths.append(p)
    # fallback to ligand.sdf if nothing else
    if not paths and Path("ligand.sdf").exists():
        paths.append(Path("ligand.sdf"))
    return paths

# ---------------------------
# Hotspot reporter (receptor-centered)
# ---------------------------
PROBE_LIBRARY = {
    "IPA":  "CC(C)O",
    "ACN":  "CC#N",
    "IMD":  "c1[nH]cnc1",
    "ACM":  "CC(=O)N",
    "PHO":  "c1ccc(cc1)O",
    "HAC":  "CC(=O)O",
    # Legacy aliases accepted elsewhere:
    "ACEA": "CC(=O)N",
    "PHOL": "c1ccc(cc1)O",
    "ACOH": "CC(=O)O",
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())

class ProbeHotspotReporter:
    """
    Write probe residue centroids (nm) every 'reportInterval' steps to CSV.
    Columns: step,resname,resid,x_nm,y_nm,z_nm

    Logs in a receptor-centered frame:
      centroid' = centroid - COM_protein(t) + COM_protein(0),
    where COM_protein uses protein Cα atoms (fallback: all protein heavy atoms).
    """
    def __init__(self, file_path, reportInterval, topology, probe_resnames=None, center_mode="protein-com"):
        self.file = open(file_path, "w")
        self.reportInterval = int(reportInterval)
        self.center_mode = center_mode
        self._have_ref_com = False
        self._ref_com = (0.0, 0.0, 0.0)

        # discover probe residue names present
        present = {(res.name or "").upper().strip() for res in topology.residues()}
        auto_names = present & PROBE_NAME_SET
        if probe_resnames:
            user_names = {s.upper().strip() for s in probe_resnames}
            self.probe_res = auto_names | (user_names & PROBE_NAME_SET)
        else:
            self.probe_res = auto_names

        # index groups to average probe centroids per residue
        groups, resid_map = [], []
        protein_resname_set = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        self._protein_anchor_indices = []
        for res in topology.residues():
            rn = (res.name or "").upper().strip()

            if rn in self.probe_res:
                idxs = [a.index for a in res.atoms()]
                if idxs:
                    groups.append(idxs)
                    resid_map.append((rn, int(res.id) if res.id is not None else 0))

            # collect protein Cα as anchors
            if rn in protein_resname_set:
                for a in res.atoms():
                    nm = (a.name or "").upper()
                    if nm == "CA":
                        self._protein_anchor_indices.append(a.index)

        # Fallback: all non-H protein atoms if no CA found
        if not self._protein_anchor_indices:
            for res in topology.residues():
                rn = (res.name or "").upper().strip()
                if rn in protein_resname_set:
                    for a in res.atoms():
                        if a.element is None or a.element.symbol == "H":
                            continue
                        self._protein_anchor_indices.append(a.index)

        self.groups = groups
        self.resid_map = resid_map

        # header + diagnostic
        self.file.write("# tracking_resnames=" + ";".join(sorted(self.probe_res)) + "\n")
        self.file.write("# center_mode=" + str(self.center_mode) + "\n")
        self.file.write("step,resname,resid,x_nm,y_nm,z_nm\n")
        self.file.flush()

    def describeNextReport(self, simulation):
        # (steps, positions, velocities, forces, energies, enforcePeriodicBox)
        return (self.reportInterval, True, False, False, False, True)

    @staticmethod
    def _com_of_indices(positions, idxs):
        sx = sy = sz = 0.0
        n = 0
        for i in idxs:
            v = positions[i]
            sx += float(v.x); sy += float(v.y); sz += float(v.z)
            n += 1
        if n == 0:
            return (0.0, 0.0, 0.0)
        return (sx/n, sy/n, sz/n)

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=False)  # list of Vec3 (nm)

        # robust step retrieval
        try:
            step = int(simulation.currentStep)  # OpenMM ≥7.5
        except Exception:
            try:
                t  = state.getTime()
                dt = simulation.integrator.getStepSize()
                step = int(round((t/dt)))
            except Exception:
                step = -1

        # receptor-centered recentering
        if self.center_mode == "protein-com" and self._protein_anchor_indices:
            com_now = self._com_of_indices(pos, self._protein_anchor_indices)
            if not self._have_ref_com:
                self._ref_com = com_now
                self._have_ref_com = True
            dx = self._ref_com[0] - com_now[0]
            dy = self._ref_com[1] - com_now[1]
            dz = self._ref_com[2] - com_now[2]
        else:
            dx = dy = dz = 0.0

        # write probe centroids
        for (idxs, (resname, resid)) in zip(self.groups, self.resid_map):
            sx = sy = sz = 0.0
            n = 0
            for i in idxs:
                v = pos[i]  # Vec3 (nm)
                sx += float(v.x); sy += float(v.y); sz += float(v.z)
                n += 1
            if n > 0:
                self.file.write(f"{step},{resname},{resid},{sx/n+dx:.5f},{sy/n+dy:.5f},{sz/n+dz:.5f}\n")
        self.file.flush()

    def __del__(self):
        try:
            self.file.close()
        except Exception:
            pass

# ---------------------------
# Packmol placements → add probes
# ---------------------------

def _default_packmol_csv(receptor_path: Path):
    root = f"{receptor_path.stem}_mixmd"
    return Path("build") / f"{root}_placements.csv"

def _openff_conf_to_angstrom_list(off_mol: Molecule):
    # Return list[(x,y,z)_Å] from OpenFF conformer without NumPy/JAX
    conf = off_mol.conformers[0]
    try:
        arr = conf.to("angstrom").m
        return [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
    except Exception:
        q = to_openmm(conf)  # Quantity(list(Vec3), angstrom)
        out = []
        for v in q.value_in_unit(u.angstrom):
            out.append((float(v[0]), float(v[1]), float(v[2])))
        return out

def _positions_for_centroid_nm(off_mol: Molecule, centroid_nm_tuple):
    coords_A = _openff_conf_to_angstrom_list(off_mol)
    n = len(coords_A)
    if n == 0:
        raise ValueError("OFF molecule has no coordinates.")

    # centroid in Å
    sx = sy = sz = 0.0
    for xA, yA, zA in coords_A:
        sx += xA; sy += yA; sz += zA
    cx = sx / n; cy = sy / n; cz = sz / n

    tx, ty, tz = (float(centroid_nm_tuple[0]),
                  float(centroid_nm_tuple[1]),
                  float(centroid_nm_tuple[2]))

    placed = []
    for xA, yA, zA in coords_A:
        xnm = (xA - cx) * 0.1 + tx
        ynm = (yA - cy) * 0.1 + ty
        znm = (zA - cz) * 0.1 + tz
        placed.append(mm.Vec3(xnm, ynm, znm))
    return [p * u.nanometer for p in placed]

def _read_centroids_csv(csv_path: Path):
    rows = []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({
                "resname": row["resname"].strip().upper(),
                "resid":   int(row["resid"]),
                "x_nm":    float(row["x_nm"]),
                "y_nm":    float(row["y_nm"]),
                "z_nm":    float(row["z_nm"]),
            })
    return rows

def _add_probes_from_packmol(modeller: Modeller,
                             forcefield: ForceField,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None = None,
                             resname_list: list[str] | None = None) -> int:
    """
    Add probes at Packmol centroids cropped to the intended MD box.
    Rename newly added residues to their probe resname (IPA/ACN/...).
    """
    if placements_csv is None:
        placements_csv = _default_packmol_csv(receptor_path)
    print(f"ℹ️  MixMD inputs: CSV={placements_csv}")
    if not placements_csv.exists():
        raise FileNotFoundError(f"Placements CSV not found: {placements_csv}")

    half = 0.5 * float(sim_box_nm)
    keep_lo = - (half - float(edge_margin_nm))
    keep_hi = + (half - float(edge_margin_nm))

    whitelist = {r.upper().strip() for r in resname_list} if resname_list else None

    rows = []
    for rec in _read_centroids_csv(placements_csv):
        resname = rec["resname"]
        if resname not in PROBE_LIBRARY:
            continue
        if whitelist and resname not in whitelist:
            continue
        x, y, z = rec["x_nm"], rec["y_nm"], rec["z_nm"]
        if (keep_lo <= x <= keep_hi) and (keep_lo <= y <= keep_hi) and (keep_lo <= z <= keep_hi):
            rows.append((resname, rec["resid"], (x, y, z)))

    if not rows:
        print("⚠️  No probe centroids inside the target box; continuing with receptor-only.")
        return 0

    used_res = sorted({res for (res, _, _) in rows})
    off_mols = []
    for res in used_res:
        off = Molecule.from_smiles(PROBE_LIBRARY[res], allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        off.name = res
        off_mols.append(off)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    top_cache = {m.name: OFFTopology.from_molecules([m]).to_openmm() for m in off_mols}
    mol_cache = {m.name: m for m in off_mols}

    def _res_list():
        return list(modeller.topology.residues())

    added = 0
    for (resname, pack_resid, cen_nm) in rows:
        top = top_cache[resname]
        off = mol_cache[resname]
        pos = _positions_for_centroid_nm(off, cen_nm)

        before_res = len(_res_list())
        modeller.add(top, pos)
        after = _res_list()

        new_res = after[before_res:]
        for res in new_res:
            try:
                res.name = resname
            except Exception:
                pass
            try:
                res.id = str(pack_resid)
            except Exception:
                pass

        added += 1

    print(f"✅ Placed {added} probe molecules (cropped to {sim_box_nm:.1f} nm box).")
    return added

# ---------------------------
# CLI
# ---------------------------

parser = argparse.ArgumentParser(description="Run OpenMM MD (RNA, protein, multi‑ligand, MixMD-from-Packmol).")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default=None, help="(Optional) single ligand SDF; otherwise auto-detect ligand_*.sdf")
parser.add_argument("--n-steps", type=int, default=1_000_000, help="Production steps (default: 1,000,000 ≈ 2 ns with 2 fs)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed")

# MixMD
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0, help="MD cubic box edge length (nm)")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACM,PHO,HAC,ACEA,PHOL,ACOH", help="Comma-separated probe residue names")
parser.add_argument("--mixmd-from-packmol", action="store_true", help="Import Packmol probes (implies --no-ligand)")
parser.add_argument("--mixmd-placements-csv", type=str, default=None, help="Override placements CSV")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15, help="Safety margin from box edge (nm)")
parser.add_argument("--hotspot-csv", type=str, default=None, help="Output CSV for probe centroids")
parser.add_argument("--hotspot-stride", type=int, default=1000, help="Stride for hotspots logging (steps)")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo
if args.mixmd_from_packmol:
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---------------------------
# Load receptor
# ---------------------------
receptor_pdb = PDBFile(args.input_receptor)

# ---------------------------
# Force field
# ---------------------------
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ---------------------------
# Build modeller & add ligands (if any)
# ---------------------------
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

ligand_paths = []
if not args.no_ligand and not args.mixmd_from_packmol:
    ligand_paths = _glob_ligands(args.input_ligand)
    if ligand_paths:
        lig_mols = []
        for lp in ligand_paths:
            m = Molecule.from_file(lp.as_posix())
            if m.n_conformers == 0:
                raise RuntimeError(f"Ligand has no conformer: {lp}")
            lig_mols.append(m)

        # Register all ligands together for SMIRNOFF
        smirnoff = SMIRNOFFTemplateGenerator(molecules=lig_mols)
        forcefield.registerTemplateGenerator(smirnoff.generator)

        # Add ligands to modeller
        for m in lig_mols:
            lig_pos = to_openmm(m.conformers[0])
            lig_top = OFFTopology.from_molecules([m]).to_openmm()
            modeller.add(lig_top, lig_pos)
        print(f"✅ Ligands merged: {', '.join(p.name for p in ligand_paths)}")
    else:
        print("ℹ️  No ligand files found (ligand.sdf or ligand_*.sdf); continuing apo.")

else:
    print("✅ No ligand: Running receptor-only (apo) setup")

# ---------------------------
# If Packmol probes are requested, add them now (before solvation)
# ---------------------------
if args.mixmd_from_packmol:
    resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    added = _add_probes_from_packmol(
        modeller=modeller,
        forcefield=forcefield,
        receptor_path=receptor_path,
        sim_box_nm=float(args.mixmd_box_size_nm),
        edge_margin_nm=float(args.mixmd_edge_margin_nm),
        placements_csv=(Path(args.mixmd_placements_csv) if args.mixmd_placements_csv else None),
        resname_list=resnames
    )
    print(f"🧪 MixMD: added {added} probe instances.")

# Hydrogens after all non‑water molecules are present
modeller.addHydrogens(forcefield)

# Save combined (pre-solvation) for reference
with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# Diagnostics before solvation
n_atoms_pre = _count_atoms(modeller.topology)
print(f"[diag] atoms pre‑solvation = {n_atoms_pre}")

# ---------------------------
# Solvate
# ---------------------------
box_edge = float(args.mixmd_box_size_nm if args.mixmd_from_packmol else 7.0)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(box_edge, box_edge, box_edge) * u.nanometer,
    ionicStrength=0.15 * u.molar,
    neutralize=True
)

with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm")

# Diagnostics after solvation
n_atoms_post = _count_atoms(modeller.topology)
print(f"[diag] atoms post‑solvation = {n_atoms_post}")

# ---------------------------
# Create System
# ---------------------------
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * u.nanometer,
    constraints=HBonds
)

seed = int(args.seed)
barostat = MonteCarloBarostat(1 * u.bar, 300 * u.kelvin, 25)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

integrator = LangevinMiddleIntegrator(300 * u.kelvin, 1 / u.picosecond, 0.002 * u.picoseconds)
integrator.setRandomNumberSeed(seed)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# --- Hotspot reporter (only if MixMD) ---
if args.mixmd_from_packmol:
    probe_resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    hotspot_csv = args.hotspot_csv or f"mixmd_hotspots{suffix}.csv"
    hrep = ProbeHotspotReporter(
        hotspot_csv,
        reportInterval=int(args.hotspot_stride),
        topology=simulation.topology,
        probe_resnames=probe_resnames
    )
    # Force an initial snapshot at step 0 so the CSV is never empty
    state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    hrep.report(simulation, state0)
    simulation.reporters.append(hrep)

# ---------------------------
# Run
# ---------------------------
print("🔹 Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter(f"minimized{suffix}.pdb", 100))

print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * u.kelvin, seed)
simulation.reporters.append(PDBReporter(f"nvt_equilibrated{suffix}.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(500)  # 1 ps

print("🔹 NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter(f"npt_equilibrated{suffix}.pdb", 500))
simulation.step(2500)  # 5 ps

print(f"🔥 Production MD: {args.n_steps:,} steps")
# Make DCD first so it always sees the same topology as final structure
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(args.n_steps)

# --- Final structures: PDB + mmCIF (to avoid serial overflow in large systems) ---
state_final = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, state_final.getPositions(), f)

# mmCIF (ChimeraX-friendly for very large atom counts)
try:
    with open(f"final_structure{suffix}.cif", "w") as f:
        PDBxFile.writeFile(simulation.topology, state_final.getPositions(), f)
    print(f"✅ Final structures saved: final_structure{suffix}.pdb and .cif")
except Exception as e:
    print(f"⚠️  mmCIF write failed ({e}); PDB written.")

# Final diagnostics
n_atoms_final = _count_atoms(simulation.topology)
print(f"[diag] atoms final(topology) = {n_atoms_final}")
print("🎉 MD complete.")
