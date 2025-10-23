# 5_md_simulation.py
# Author: Iori Mochizuki + collaborator patch
# Updated: 2025-10-23
# Description: OpenMM 2 ns MD (protein/RNA/ligand) with optional MixMD-from-Packmol import.
# - MixMD (Packmol) implies apo; probes added before solvation and logged every stride to CSV
# - Multi-ligand: auto-detect ligand*.sdf and add all (or fall back to --input-ligand)
# - Writes both PDB and mmCIF for large systems to avoid PDB serial overflow issues

import argparse, csv, glob
from pathlib import Path
from sys import stdout

import openmm as mm
from openmm.app import *
from openmm import MonteCarloBarostat, LangevinMiddleIntegrator
from openmm.unit import *
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from openff.toolkit.topology import Molecule, Topology as OFFTopology
from openff.units.openmm import to_openmm

# ---------------------------
# Hotspot reporter (manual reporting; receptor-centered)
# ---------------------------
class ProbeHotspotReporter:
    """
    Write probe residue centroids (nm) to CSV.
    Columns: step,resname,resid,x_nm,y_nm,z_nm

    We log in a receptor-centered frame:
      centroid' = centroid - COM_protein(t) + COM_protein(0),
    where COM_protein uses protein Cα atoms (fallback: protein heavy atoms).
    """
    def __init__(self, file_path, topology, probe_resnames=None, center_mode="protein-com"):
        self.file = open(file_path, "w")
        self.center_mode = center_mode
        self._have_ref_com = False
        self._ref_com = (0.0, 0.0, 0.0)

        # Discover probe residue names present
        present = {(res.name or "").upper().strip() for res in topology.residues()}
        auto_names = present & PROBE_NAME_SET
        if probe_resnames:
            user_names = {s.upper().strip() for s in probe_resnames}
            self.probe_res = auto_names | (user_names & PROBE_NAME_SET)
        else:
            self.probe_res = auto_names

        # Build index groups to average probe centroids per residue
        groups, resid_map = [], []
        protein_resname_set = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        self._protein_anchor_indices = []
        for res in topology.residues():
            rn = (res.name or "").upper().strip()

            # Track probe residues
            if rn in self.probe_res:
                idxs = [a.index for a in res.atoms()]
                if idxs:
                    groups.append(idxs)
                    resid_map.append((rn, int(res.id) if res.id is not None else 0))

            # Collect anchors (Cα preferred)
            if rn in protein_resname_set:
                for a in res.atoms():
                    nm = (a.name or "").upper()
                    if nm == "CA":
                        self._protein_anchor_indices.append(a.index)

        # Fallback: heavy atoms if no Cα
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

    def write(self, simulation, step, state):
        """Manual logging entry point (call whenever you want a row)."""
        pos = state.getPositions(asNumpy=False)  # list of Vec3 (nm)

        # receptor-centered recentering
        if self._protein_anchor_indices:
            com_now = self._com_of_indices(pos, self._protein_anchor_indices)
            if not self._have_ref_com:
                self._ref_com = com_now
                self._have_ref_com = True
            dx = self._ref_com[0] - com_now[0]
            dy = self._ref_com[1] - com_now[1]
            dz = self._ref_com[2] - com_now[2]
        else:
            dx = dy = dz = 0.0

        for (idxs, (resname, resid)) in zip(self.groups, self.resid_map):
            sx = sy = sz = 0.0
            n = 0
            for i in idxs:
                v = pos[i]
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


# --- MixMD-from-Packmol helpers (pure OpenMM, no NumPy) -----------------------
# Canonical 3-letter names used by 3c/Packmol, plus 4-letter aliases so 5 accepts both.
PROBE_LIBRARY = {
    "IPA":  "CC(C)O",           # isopropanol
    "ACN":  "CC#N",             # acetonitrile
    "IMD":  "c1[nH]cnc1",       # imidazole (neutral, N1–H)
    "ACM":  "CC(=O)N",          # acetamide
    "PHO":  "c1ccc(cc1)O",      # phenol
    "HAC":  "CC(=O)O",          # acetic acid
    # legacy aliases accepted on input:
    "ACEA": "CC(=O)N", "PHOL": "c1ccc(cc1)O", "ACOH": "CC(=O)O",
}
PROBE_NAME_SET = set(PROBE_LIBRARY.keys())

def _default_packmol_paths(receptor_path: Path):
    root = f"{receptor_path.stem}_mixmd"
    placements  = Path("build") / f"{root}_placements.csv"
    return placements

def _openff_conf_to_angstrom_list(off_mol):
    """Return a list of (x,y,z) floats in Å from an OpenFF Molecule conformer."""
    conf = off_mol.conformers[0]  # OpenFF quantity w/ units
    try:
        arr = conf.to("angstrom").m  # N×3 numeric
        return [(float(r[0]), float(r[1]), float(r[2])) for r in arr]
    except Exception:
        q = to_openmm(conf)  # Quantity(list(Vec3), angstrom)
        out = []
        for v in q.value_in_unit(angstrom):
            out.append((float(v[0]), float(v[1]), float(v[2])))
        return out

def _positions_for_centroid_nm(off_mol, centroid_nm_tuple):
    """Build OpenMM positions (list[Vec3]*nanometer) centered on centroid_nm_tuple (nm)."""
    coords_A = _openff_conf_to_angstrom_list(off_mol)
    n = len(coords_A)
    if n == 0:
        raise ValueError("OFF molecule has no coordinates.")

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
    return [p * nanometer for p in placed]

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

def _add_probes_from_packmol(modeller,
                             forcefield,
                             receptor_path: Path,
                             sim_box_nm: float,
                             edge_margin_nm: float,
                             placements_csv: Path | None = None,
                             resname_list: list[str] | None = None) -> int:
    """Add probes at Packmol centroids cropped to the intended MD box; rename residues."""
    if placements_csv is None:
        placements_csv = _default_packmol_paths(receptor_path)

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
        smi = PROBE_LIBRARY[res]
        off = Molecule.from_smiles(smi, allow_undefined_stereo=True)
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
parser = argparse.ArgumentParser(description="Run OpenMM MD with optional MixMD-from-Packmol and multi-ligand support.")
parser.add_argument("--rna", action="store_true", help="Enable RNA forcefield and logic")
parser.add_argument("--no-ligand", action="store_true", help="Run without ligand (receptor only)")
parser.add_argument("--input-receptor", type=str, default="receptor_cleaned.pdb", help="Input receptor PDB")
parser.add_argument("--input-ligand", type=str, default="ligand.sdf", help="Single-ligand fallback if no ligand*.sdf detected")
parser.add_argument("--n-steps", type=int, default=1000000, help="Production MD steps (default: 1,000,000 ~ 2 ns at 2 fs)")
parser.add_argument("--seed", type=int, default=13579, help="Random seed for barostat/integrator/velocities")

# Multi‑ligand discovery
parser.add_argument("--ligand-glob", type=str, default="ligand*.sdf",
                    help="Glob pattern to discover multiple ligands (default: ligand*.sdf).")

# MixMD-from-Packmol options
parser.add_argument("--mixmd-box-size-nm", type=float, default=7.0,
                    help="MD cubic box edge length; probes outside this cube are dropped.")
parser.add_argument("--mixmd-resnames", type=str, default="IPA,ACN,IMD,ACM,PHO,HAC,ACEA,PHOL,ACOH",
                    help="Comma-separated probe residue names to import.")
parser.add_argument("--mixmd-from-packmol", action="store_true",
    help="Read placements CSV (build/<receptor>_mixmd_placements.csv) and place probes before solvation (implies --no-ligand).")
parser.add_argument("--mixmd-placements-csv", type=str, default=None,
    help="Override path to placements CSV (default: build/<receptor>_mixmd_placements.csv).")
parser.add_argument("--mixmd-edge-margin-nm", type=float, default=0.15,
    help="Safety margin from box edge when filtering centroids.")
parser.add_argument("--hotspot-csv", type=str, default=None,
    help="If set with --mixmd-from-packmol, write probe centroids every --hotspot-stride to this CSV (default: mixmd_hotspots_no_ligand.csv).")
parser.add_argument("--hotspot-stride", type=int, default=1000,
    help="Stride (steps) for hotspot centroid logging.")

# Large systems output
parser.add_argument("--write-cif", action="store_true",
    help="Also write mmCIF versions of solvated and final structures (recommended for large boxes).")

args = parser.parse_args()
receptor_path = Path(args.input_receptor)

# MixMD implies apo
if getattr(args, "mixmd_from_packmol", False):
    args.no_ligand = True

suffix = "_no_ligand" if args.no_ligand else ""

# ---------------------------
# Load receptor
# ---------------------------
receptor_pdb = PDBFile(args.input_receptor)

# ---------------------------
# Force fields
# ---------------------------
if args.rna:
    print("🧬 [RNA MODE] Using amber14/RNA.OL3 + TIP3P-FB")
    forcefield = ForceField("amber14/RNA.OL3.xml", "amber14/tip3pfb.xml")
else:
    print("🧬 [Protein MODE] Using amber14-all + TIP3P-FB")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

# ---------------------------
# Build modeller (ligands or apo)
# ---------------------------
modeller = Modeller(receptor_pdb.topology, receptor_pdb.positions)

lig_files = []
if not args.no_ligand and not args.mixmd_from_packmol:
    lig_files = sorted(glob.glob(args.ligand_glob))
    if not lig_files and Path(args.input_ligand).exists():
        lig_files = [args.input_ligand]

if lig_files and not args.no_ligand:
    print(f"🔗 Detected {len(lig_files)} ligand SDF file(s): {lig_files}")
    lig_mols = []
    for path in lig_files:
        m = Molecule.from_file(path)
        if m.n_conformers == 0:
            raise ValueError(f"Ligand has no conformer: {path}")
        m.name = Path(path).stem[:3].upper()  # readable short tag
        lig_mols.append(m)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=lig_mols)
    forcefield.registerTemplateGenerator(smirnoff.generator)
    for m in lig_mols:
        lig_top = OFFTopology.from_molecules([m]).to_openmm()
        lig_pos = to_openmm(m.conformers[0])
        modeller.add(lig_top, lig_pos)
    print("✅ Ligand(s) merged")
else:
    if args.no_ligand:
        print("✅ No ligand: Running receptor-only (apo) setup")
    else:
        print("⚠️  No ligand SDF detected; continuing without ligands (apo).")
        args.no_ligand = True
        suffix = "_no_ligand"

# If MixMD probes are requested, add them now (before solvation)
if args.mixmd_from_packmol:
    resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    added = _add_probes_from_packmol(
        modeller=modeller,
        forcefield=forcefield,
        receptor_path=receptor_path,
        sim_box_nm=float(args.mixmd_box_size_n m),  # <-- ensure no stray space in your editor
        edge_margin_nm=float(args.mixmd_edge_margin_nm),
        placements_csv=(Path(args.mixmd_placements_csv) if args.mixmd_placements_csv else None),
        resname_list=resnames
    )
    print(f"🧪 MixMD: added {added} probe instances.")

# Hydrogens after all non‑water molecules are present
modeller.addHydrogens(forcefield)

with open(f"combined_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ System ready for solvation ({'apo' if args.no_ligand else 'holo'})")

# ---------------------------
# Solvate in the requested MD box
# ---------------------------
box_edge = float(args.mixmd_box_size_nm if args.mixmd_from_packmol else 7.0)
modeller.addSolvent(
    forcefield,
    model="tip3p",
    boxSize=(box_edge, box_edge, box_edge) * nanometer,
    ionicStrength=0.15 * molar,
    neutralize=True
)
with open(f"solvated_receptor_ligand{suffix}.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
if args.write_cif:
    # mmCIF avoids PDB serial overflow for large systems
    with open(f"solvated_receptor_ligand{suffix}.cif", "w") as f:
        PDBxFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✅ Solvated system ready. Box = {box_edge:.2f} nm")

# ---------------------------
# Create System
# ---------------------------
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * nanometer,
    constraints=HBonds
)

seed = int(args.seed)
barostat = MonteCarloBarostat(1 * bar, 300 * kelvin, 25)
barostat.setRandomNumberSeed(seed)
system.addForce(barostat)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
integrator.setRandomNumberSeed(seed)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# ---------------------------
# Hotspot logging (MixMD only): manual loop
# ---------------------------
hrep = None
if args.mixmd_from_packmol:
    probe_resnames = [s.strip().upper() for s in args.mixmd_resnames.split(",") if s.strip()]
    hotspot_csv = args.hotspot_csv or f"mixmd_hotspots{suffix}.csv"
    hrep = ProbeHotspotReporter(
        hotspot_csv,
        topology=simulation.topology,
        probe_resnames=probe_resnames
    )
    # initial snapshot at step 0
    state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    hrep.write(simulation, step=0, state=state0)

def _manual_run(simulation, total_steps: int, stride: int):
    """Run `total_steps` in chunks of `stride`, manually logging hotspots if enabled."""
    remaining = int(total_steps)
    while remaining > 0:
        n = min(int(stride), remaining)
        simulation.step(n)
        remaining -= n
        if hrep is not None:
            # robust step retrieval
            try:
                step = int(simulation.currentStep)
            except Exception:
                st = simulation.context.getState(getTime=True)
                dt = simulation.integrator.getStepSize()
                step = int(round((st.getTime()/dt)))
            state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
            hrep.write(simulation, step=step, state=state)

# ---------------------------
# Run
# ---------------------------
print("🔹 Energy Minimization...")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter(f"minimized{suffix}.pdb", 200))

print("🔹 NVT Equilibration (1 ps)...")
simulation.context.setVelocitiesToTemperature(300 * kelvin, seed)
simulation.reporters.append(PDBReporter(f"nvt_equilibrated{suffix}.pdb", 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
# Run NVT with manual hotspot logging if enabled
_manual_run(simulation, total_steps=500, stride=max(1, int(args.hotspot_stride)))

print("🔹 NPT Equilibration (5 ps)...")
simulation.reporters.append(PDBReporter(f"npt_equilibrated{suffix}.pdb", 500))
_manual_run(simulation, total_steps=2500, stride=max(1, int(args.hotspot_stride)))

print(f"🔥 Production MD: {args.n_steps:,} steps")
simulation.reporters.append(DCDReporter(f"production_md{suffix}.dcd", 1000))
simulation.reporters.append(StateDataReporter(f"production_md{suffix}.log", 1000, step=True, potentialEnergy=True, temperature=True))
_manual_run(simulation, total_steps=int(args.n_steps), stride=max(1, int(args.hotspot_stride)))

# Final coordinates
state_final = simulation.context.getState(getPositions=True)
with open(f"final_structure{suffix}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, state_final.getPositions(), f)
if args.write_cif:
    with open(f"final_structure{suffix}.cif", "w") as f:
        PDBxFile.writeFile(simulation.topology, state_final.getPositions(), f)

print(f"🎉 MD complete → final_structure{suffix}.pdb saved.")
if args.write_cif:
    print(f"📦 Also wrote mmCIF → final_structure{suffix}.cif (recommended for large boxes)")
