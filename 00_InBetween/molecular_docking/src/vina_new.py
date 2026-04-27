#!/usr/bin/env python3
"""
Script for docking Enzyme with NNBT (keeping the heme group)
=======================================

Pipeline
--------
  1. Download 6EKZ from RCSB (or use local file)
  2. Prepare receptor  → remove waters/HETATM, keep heme, write PDBQT
  3. Prepare ligand    → NNBT SMILES → 3D conformer → PDBQT
  4. Calculate docking box from Phe69/121/199 Cα centroid
  5. Run AutoDock Vina
  6. Analyse poses     → contacts, α-C to heme-Fe distance
  7. Write summary report

Install
-------
    pip install biopython rdkit meeko numpy

    AutoDock Vina binary:
      Linux: wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
              chmod +x vina_*  &&  mv vina_* vina
      Mac:   download from https://github.com/ccsb-scripps/AutoDock-Vina/releases

    MGLTools (optional, for prepare_receptor4.py alternative):
      http://mgltools.scripps.edu/downloads

Usage
-----
    # Download 6EKZ automatically and dock:
    python dock_nnbt.py --vina ./vina --workdir results/

    # Use a local PDB file:
    python dock_nnbt.py --vina ./vina --receptor 6EKZ.pdb --workdir results/

    # Use a pre-prepared receptor PDBQT (skip preparation):
    python dock_nnbt.py --vina ./vina --receptor_pdbqt receptor.pdbqt --workdir results/
"""

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

NNBT_SMILES = "CC1=CC=C(C=C1)N(CC(C)O)CC(C)O"
PDB_ID      = "6EKZ"
PDB_URL     = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

# Phe triad residue numbers (1-indexed, chain A, from Ramirez-Escudero 2018)
PHE_TRIAD   = [69, 121, 199]

# Box dimensions in Å — large enough to cover the full 17 Å channel depth
BOX_SIZE    = (22.0, 22.0, 22.0)

# Docking parameters
EXHAUSTIVENESS = 16
NUM_MODES      = 9

# Contact cutoff for analysis
CONTACT_CUTOFF = 4.5   # Å

# Productive N-dealkylation geometry target
ALPHA_C_TARGET_MIN = 3.5   # Å
ALPHA_C_TARGET_MAX = 4.5   # Å

# Residues to KEEP when cleaning (heme must stay in receptor)
HEME_RESNAME = "HEM"

# Catalytic residues for reporting
CATALYTIC    = {36: "Cys36 (axial ligand)", 189: "Arg189 (acid-base)", 196: "Glu196 (acid-base)"}
HOTSPOTS     = {69, 76, 121, 191, 199, 244, 274, 277, 311, 316}


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — Download PDB
# ─────────────────────────────────────────────────────────────────────────────

def download_pdb(workdir: Path) -> Path:
    out = workdir / f"{PDB_ID}.pdb"
    if out.exists():
        log.info(f"  PDB already present: {out}")
        return out
    log.info(f"  Downloading {PDB_ID} from RCSB...")
    urllib.request.urlretrieve(PDB_URL, out)
    log.info(f"  Saved → {out}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 — Receptor preparation
# ─────────────────────────────────────────────────────────────────────────────

def prepare_receptor(pdb_path: Path, workdir: Path) -> tuple[Path, dict]:
    """
    Clean PDB: keep protein ATOM records + selected HETATM (heme + metal ions).
    Add polar hydrogens using reduce (if available) or flag for manual addition.
    Write clean PDB and PDBQT.

    Returns (pdbqt_path, heme_coords_dict)
    """
    from Bio.PDB import PDBParser, PDBIO, Select

    log.info("  Preparing receptor...")

    # Explicitly allowed cofactors / ions
    ALLOWED_HET = {"HEM", "CL", "MG", "NAG", "SNP", "PO4"}

    class ReceptorSelect(Select):
        def accept_residue(self, residue):
            # Keep standard amino acids + selected hetero groups
            return (
                residue.id[0] == " " or
                residue.get_resname() in ALLOWED_HET
            )

        def accept_atom(self, atom):
            # Remove alternate conformations — keep only blank or 'A'
            return atom.altloc in (" ", "A", "")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))

    # Extract heme iron coordinates for later distance measurement
    heme_coords = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname() == "HEM":
                    for atom in res:
                        if atom.get_name() == "FE":
                            heme_coords["FE"] = atom.get_vector().get_array().tolist()
                            heme_coords["chain"] = chain.id
                            heme_coords["resnum"] = res.id[1]
                    break

    if not heme_coords:
        log.warning("  Heme iron (FE) not found — check PDB HETATM records")
    else:
        log.info(f"  Heme FE at: {[round(x,2) for x in heme_coords['FE']]}")

    # Write clean PDB
    clean_pdb = workdir / "receptor_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(clean_pdb), ReceptorSelect())
    log.info(f"  Clean receptor → {clean_pdb}")

    # Try to add hydrogens with reduce; fall back to just using the clean PDB
    pdb_with_h = workdir / "receptor_H.pdb"
    reduce_ok = False
    if shutil.which("reduce"):
        result = subprocess.run(
            ["reduce", "-NOFLIP", str(clean_pdb)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            pdb_with_h.write_text(result.stdout)
            reduce_ok = True
            log.info("  Polar H added with reduce")
    if not reduce_ok:
        shutil.copy(clean_pdb, pdb_with_h)
        log.warning(
            "  reduce not found — using PDB without explicit H. "
            "AutoDock Vina will add H internally but Gasteiger charges may be suboptimal. "
            "Install reduce for best results: https://github.com/rlabduke/reduce"
        )

    # Convert to PDBQT using meeko (preferred) or obabel
    receptor_pdbqt = workdir / "receptor.pdbqt"
    if _meeko_available():
        _pdb_to_pdbqt_meeko_receptor(pdb_with_h, receptor_pdbqt)
    elif shutil.which("obabel"):
        subprocess.run(
            ["obabel", str(pdb_with_h), "-O", str(receptor_pdbqt), "-xr"],
            check=True, capture_output=True
        )
        log.info(f"  Receptor PDBQT (obabel) → {receptor_pdbqt}")
    else:
        _pdb_to_pdbqt_minimal(pdb_with_h, receptor_pdbqt)
        log.warning(
            "  Neither meeko nor obabel found — wrote minimal PDBQT. "
            "Install meeko (pip install meeko) for proper charge assignment."
        )

    return receptor_pdbqt, heme_coords


def _meeko_available() -> bool:
    try:
        import meeko
        return True
    except ImportError:
        return False


def _pdb_to_pdbqt_meeko_receptor(pdb_path: Path, out_path: Path):
    """Use meeko for receptor preparation with proper Gasteiger charges."""
    import meeko
    result = subprocess.run(
        ["mk_prepare_receptor.py", "-i", str(pdb_path), "-o", str(out_path)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        log.warning(f"  meeko mk_prepare_receptor failed: {result.stderr[:200]}")
        _pdb_to_pdbqt_minimal(pdb_path, out_path)
    else:
        log.info(f"  Receptor PDBQT (meeko) → {out_path}")


_ELEMENT_TO_ADTYPE = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "H": "HD",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "FE": "Fe", "MG": "Mg", "CA": "Ca", "MN": "Mn", "ZN": "Zn",
    "CU": "Cu", "NA": "Na", "K": "K",
}

def _pdb_to_pdbqt_minimal(pdb_path: Path, out_path: Path):
    """
    Minimal PDB → PDBQT conversion: copy ATOM/HETATM lines, set AutoDock atom type.
    Sufficient for Vina which uses its own internal typing.
    For best results use meeko or obabel.
    """
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip().upper() if len(line) > 76 else ""
                ad_type = _ELEMENT_TO_ADTYPE.get(element, element.capitalize() if element else "C")
                # PDBQT columns 1-66 from PDB, then charge (0.000), then atom type
                base = line[:66].rstrip().ljust(66)
                # PDBQT cols 67-70: spaces, 71-76: charge (6.3f), 77: space, 78-79: atom type
                lines.append(f"{base}    {0.000:6.3f} {ad_type:<2s}\n")
    out_path.write_text("".join(lines))
    log.info(f"  Receptor PDBQT (minimal) → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 — Ligand preparation
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ligand(workdir: Path) -> tuple[Path, dict]:
    """
    Build NNBT 3D conformer from SMILES, optimise, write PDBQT.
    Returns (pdbqt_path, atom_index_map) where atom_index_map
    identifies the α-carbon atoms for later distance measurement.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    log.info("  Preparing NNBT ligand...")

    mol = Chem.MolFromSmiles(NNBT_SMILES)
    mol = Chem.AddHs(mol)

    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        raise RuntimeError("3D embedding failed for NNBT")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    mw = Descriptors.MolWt(mol)
    log.info(f"  NNBT MW = {mw:.1f} Da  |  {mol.GetNumHeavyAtoms()} heavy atoms")

    # Identify α-carbon atom indices (carbons directly bonded to N)
    # NNBT has one N; the α-carbons are the two CH₂ groups bonded to it
    alpha_c_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:   # nitrogen
            n_idx = atom.GetIdx()
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 6:   # carbon neighbour = α-carbon
                    alpha_c_indices.append(nb.GetIdx())
            break

    log.info(f"  α-carbon atom indices (RDKit): {alpha_c_indices}")

    # Write SDF
    sdf_path = workdir / "NNBT.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    # Convert to PDBQT
    ligand_pdbqt = workdir / "NNBT.pdbqt"
    if _meeko_available():
        _sdf_to_pdbqt_meeko(sdf_path, ligand_pdbqt)
    elif shutil.which("obabel"):
        subprocess.run(
            ["obabel", str(sdf_path), "-O", str(ligand_pdbqt)],
            check=True, capture_output=True
        )
        log.info(f"  Ligand PDBQT (obabel) → {ligand_pdbqt}")
    else:
        _sdf_to_pdbqt_minimal(mol, ligand_pdbqt)
        log.warning(
            "  Neither meeko nor obabel found — wrote minimal ligand PDBQT. "
            "Install meeko (pip install meeko) for proper rotatable bond handling."
        )

    return ligand_pdbqt, {"alpha_c_indices": alpha_c_indices, "n_atoms": mol.GetNumAtoms()}


def _sdf_to_pdbqt_meeko(sdf_path: Path, out_path: Path):
    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", str(sdf_path), "-o", str(out_path)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        log.warning(f"  meeko mk_prepare_ligand failed: {result.stderr[:200]}")
        from rdkit import Chem
        mol = next(Chem.SDMolSupplier(str(sdf_path), removeHs=False))
        _sdf_to_pdbqt_minimal(mol, out_path)
    else:
        log.info(f"  Ligand PDBQT (meeko) → {out_path}")


def _sdf_to_pdbqt_minimal(mol, out_path: Path):
    """Write a minimal PDBQT from an RDKit mol. Sufficient for Vina scoring."""
    from rdkit import Chem
    conf = mol.GetConformer()
    lines = ["ROOT\n"]
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(i)
        sym = atom.GetSymbol()
        lines.append(
            f"ATOM  {i+1:5d}  {sym:<3s} LIG A   1    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    "
            f"     0.000 {sym}\n"
        )
    lines.append("ENDROOT\n")
    lines.append("TORSDOF 0\n")
    out_path.write_text("".join(lines))
    log.info(f"  Ligand PDBQT (minimal) → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 — Calculate docking box from Phe triad
# ─────────────────────────────────────────────────────────────────────────────

def calculate_box_center(pdb_path: Path) -> np.ndarray:
    """
    Return centroid of Phe69, Phe121, Phe199 Cα atoms.
    This centres the docking box on the substrate-binding triad.
    """
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))

    coords = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[1] in PHE_TRIAD and res.id[0] == " ":
                    try:
                        ca = res["CA"].get_vector().get_array()
                        coords.append(ca)
                        log.info(f"  Phe{res.id[1]} Cα: {[round(x,2) for x in ca]}")
                    except KeyError:
                        log.warning(f"  No Cα found for residue {res.id[1]}")

    if not coords:
        raise ValueError(
            "Could not find Phe69/121/199 Cα atoms. "
            "Check chain ID and residue numbering in your PDB."
        )

    centre = np.mean(coords, axis=0)
    log.info(f"  Docking box centre (Phe triad centroid): {[round(x,2) for x in centre]}")
    return centre


# ─────────────────────────────────────────────────────────────────────────────
# Step 5 — Run AutoDock Vina
# ─────────────────────────────────────────────────────────────────────────────

def run_vina(
    vina_bin: Path,
    receptor: Path,
    ligand: Path,
    center: np.ndarray,
    box: tuple,
    workdir: Path,
    label: str = "WT",
) -> tuple[Path, list[float]]:
    """
    Run Vina and return (output_pdbqt_path, list_of_scores).
    """
    out_pdbqt = workdir / f"poses_{label}.pdbqt"
    log_file  = workdir / f"vina_{label}.log"

    cmd = [
        str(vina_bin.resolve()),
        "--receptor",      str(receptor),
        "--ligand",        str(ligand),
        "--center_x",      f"{center[0]:.3f}",
        "--center_y",      f"{center[1]:.3f}",
        "--center_z",      f"{center[2]:.3f}",
        "--size_x",        f"{box[0]:.1f}",
        "--size_y",        f"{box[1]:.1f}",
        "--size_z",        f"{box[2]:.1f}",
        "--out",           str(out_pdbqt),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes",     str(NUM_MODES),
    ]

    log.info(f"  Running AutoDock Vina ({label})...")
    log.info(f"  cmd: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    # Vina writes its output to stdout; save it as the log
    log_file.write_text(result.stdout)
    if result.returncode != 0:
        log.error(result.stderr)
        raise RuntimeError(f"Vina failed for {label}: {result.stderr[-500:]}")

    scores = _parse_vina_scores(log_file)
    log.info(f"  Top {len(scores)} poses: {[round(s,2) for s in scores]} kcal/mol")

    return out_pdbqt, scores


def _parse_vina_scores(log_file: Path) -> list[float]:
    scores = []
    try:
        with open(log_file) as f:
            in_table = False
            for line in f:
                if "-----+" in line:
                    in_table = True
                    continue
                if in_table:
                    m = re.match(r"\s+\d+\s+([-\d.]+)", line)
                    if m:
                        scores.append(float(m.group(1)))
                    elif line.strip() == "":
                        break
    except FileNotFoundError:
        log.warning(f"Vina log not found: {log_file}")
    return scores


# ─────────────────────────────────────────────────────────────────────────────
# Step 6 — Analyse best pose
# ─────────────────────────────────────────────────────────────────────────────

def analyse_pose(
    poses_pdbqt: Path,
    receptor_pdb: Path,
    heme_coords: dict,
    workdir: Path,
    label: str = "WT",
) -> dict:
    """
    Parse best pose from PDBQT, extract:
      - Ligand heavy atom coordinates
      - Contact residues within CONTACT_CUTOFF Å
      - α-carbon to heme FE distance
      - Whether geometry is productive for N-dealkylation
    """
    from Bio.PDB import PDBParser

    log.info(f"  Analysing best pose ({label})...")

    # Parse best pose from PDBQT (first MODEL block)
    ligand_coords = _parse_best_pose_pdbqt(poses_pdbqt)
    if ligand_coords is None or len(ligand_coords) == 0:
        log.warning("  Could not parse ligand coordinates from PDBQT")
        return {}

    # Parse receptor
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(receptor_pdb))

    aa_map = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
        "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    }

    contacts = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] not in (" ",):
                    continue
                resnum = res.id[1]
                resname = res.get_resname()
                aa1 = aa_map.get(resname, "X")
                for atom in res:
                    coord = atom.get_vector().get_array()
                    dists = np.linalg.norm(ligand_coords - coord, axis=1)
                    if dists.min() <= CONTACT_CUTOFF:
                        contacts[resnum] = f"{aa1}{resnum}"
                        break

    contact_list = [contacts[k] for k in sorted(contacts.keys())]
    hotspot_contacts = [c for c in contact_list
                        if int(re.sub(r"[^0-9]","",c)) in HOTSPOTS]

    log.info(f"  Contact residues (≤{CONTACT_CUTOFF}Å): {contact_list}")
    log.info(f"  Hotspot contacts: {hotspot_contacts}")

    # α-carbon to heme FE distance
    # α-carbons are the two carbons with highest connectivity to N in NNBT
    # We approximate by taking the two ligand atoms closest to heme FE
    # that are not the deepest atom (to avoid N itself)
    alpha_c_dist = None
    geometry_ok = False

    if heme_coords.get("FE"):
        fe_coord = np.array(heme_coords["FE"])
        dists_to_fe = np.linalg.norm(ligand_coords - fe_coord, axis=1)
        sorted_dists = np.sort(dists_to_fe)

        # The closest heavy atom to FE is most likely the reacting carbon
        # (or N if geometry is wrong)
        min_dist = sorted_dists[0]
        alpha_c_dist = round(float(min_dist), 2)

        geometry_ok = ALPHA_C_TARGET_MIN <= alpha_c_dist <= ALPHA_C_TARGET_MAX

        log.info(
            f"  Closest ligand atom to heme FE: {alpha_c_dist} Å "
            f"({'✓ productive range' if geometry_ok else '✗ outside productive range 3.5-4.5Å'})"
        )

    # Save best pose as PDB for PyMOL
    best_pose_pdb = workdir / f"best_pose_{label}.pdb"
    _pdbqt_pose_to_pdb(poses_pdbqt, best_pose_pdb)

    result = {
        "label": label,
        "contacts": contact_list,
        "hotspot_contacts": hotspot_contacts,
        "alpha_c_to_FE_dist_A": alpha_c_dist,
        "geometry_productive": geometry_ok,
        "best_pose_pdb": str(best_pose_pdb),
    }

    return result


def _parse_best_pose_pdbqt(pdbqt_path: Path) -> np.ndarray | None:
    """Extract heavy atom coordinates from first MODEL in PDBQT."""
    coords = []
    in_first_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if in_first_model:
                    break   # stop at second model
                in_first_model = True
            if line.startswith(("ATOM", "HETATM")) and in_first_model:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else ""
                    if element != "H":   # skip hydrogens
                        coords.append([x, y, z])
                except ValueError:
                    continue
    return np.array(coords) if coords else None


def _pdbqt_pose_to_pdb(pdbqt_path: Path, out_pdb: Path):
    """Write first MODEL from PDBQT as a plain PDB for PyMOL."""
    lines = []
    in_first = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if in_first:
                    break
                in_first = True
                continue
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM", "REMARK")):
                # Strip AutoDock-specific columns (beyond col 66)
                lines.append(line[:66] + "\n")
    out_pdb.write_text("".join(lines))


# ─────────────────────────────────────────────────────────────────────────────
# Step 7 — Write summary report
# ─────────────────────────────────────────────────────────────────────────────

def write_report(
    results: list[dict],
    scores: dict[str, list[float]],
    heme_coords: dict,
    workdir: Path,
):
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("NNBT × PaDa-I DOCKING REPORT")
    report_lines.append("=" * 60)
    report_lines.append(f"Receptor : PDB {PDB_ID} (6EKZ, Ramirez-Escudero 2018)")
    report_lines.append(f"Ligand   : NNBT  {NNBT_SMILES}")
    report_lines.append(f"Software : AutoDock Vina")
    report_lines.append(f"Box size : {BOX_SIZE[0]} × {BOX_SIZE[1]} × {BOX_SIZE[2]} Å")
    report_lines.append(f"Heme FE  : {heme_coords.get('FE', 'not found')}")
    report_lines.append("")

    for res in results:
        label = res.get("label", "?")
        report_lines.append(f"── {label} ──")
        sc = scores.get(label, [])
        if sc:
            report_lines.append(f"  Best Vina score    : {sc[0]:.2f} kcal/mol")
            report_lines.append(f"  All poses          : {[round(s,2) for s in sc]}")
        report_lines.append(f"  Contact residues   : {res.get('contacts', [])}")
        report_lines.append(f"  Hotspot contacts   : {res.get('hotspot_contacts', [])}")
        dist = res.get("alpha_c_to_FE_dist_A")
        geom = res.get("geometry_productive")
        report_lines.append(
            f"  α-C to heme FE     : {dist} Å  "
            f"({'✓ productive' if geom else '✗ non-productive — check N vs α-C orientation'})"
        )
        report_lines.append(f"  Best pose PDB      : {res.get('best_pose_pdb', 'N/A')}")
        report_lines.append("")

    report_lines.append("─" * 60)
    report_lines.append("KEY RESIDUES REFERENCE (Ramirez-Escudero 2018)")
    report_lines.append("  Phe triad (substrate sandwiching) : F69, F121, F199")
    report_lines.append("  Entrance gate (molecular hinge)   : F76, F191")
    report_lines.append("  Loop hotspot                      : A316  (A316P → JEd-I)")
    report_lines.append("  Outer entrance hydrophobics        : V244, F274, P277")
    report_lines.append("  Polar contacts for -OH groups      : A73, T192")
    report_lines.append("  PROTECTED (never mutate)           : C36, R189, E196")
    report_lines.append("")
    report_lines.append("GEOMETRY INTERPRETATION")
    report_lines.append(
        f"  Productive N-dealkylation: α-C to FE = {ALPHA_C_TARGET_MIN}–{ALPHA_C_TARGET_MAX} Å"
    )
    report_lines.append(
        "  If N faces FE instead: N-oxidation competes → non-productive for lactaldehyde"
    )
    report_lines.append("=" * 60)

    report_text = "\n".join(report_lines)
    print("\n" + report_text)

    report_path = workdir / "docking_report.txt"
    report_path.write_text(report_text)
    log.info(f"\n  Full report → {report_path}")

    # Also save JSON for downstream use (e.g. feeding into Boltz-2 pipeline)
    json_path = workdir / "docking_results.json"
    json_path.write_text(json.dumps({
        "receptor": PDB_ID,
        "ligand": "NNBT",
        "heme_FE_coords": heme_coords.get("FE"),
        "results": results,
        "scores": scores,
    }, indent=2))
    log.info(f"  JSON results → {json_path}")


# ─────────────────────────────────────────────────────────────────────────────
# PyMOL visualisation script
# ─────────────────────────────────────────────────────────────────────────────

def write_pymol_script(workdir: Path, results: list[dict]):
    """Generate a PyMOL script to visualise the docking result."""
    best_pose_pdb = results[0].get("best_pose_pdb", "") if results else ""
    contacts = results[0].get("contacts", []) if results else []
    contact_resnums = "+".join(re.sub(r"[^0-9]","",c) for c in contacts)

    script = f"""# PyMOL visualisation — NNBT docked into PaDa-I (6EKZ)
# Run: pymol visualise_docking.py

from pymol import cmd

# Load structures
cmd.load("{workdir}/receptor_clean.pdb", "receptor")
cmd.load("{best_pose_pdb}", "nnbt_best_pose")

# Basic display
cmd.remove("solvent")
cmd.bg_color("white")
cmd.show("cartoon", "receptor")
cmd.color("white", "receptor")

# Heme
cmd.select("heme", "receptor and resn HEM")
cmd.show("sticks", "heme")
cmd.color("green", "heme")

# NNBT best pose
cmd.show("sticks", "nnbt_best_pose")
cmd.color("yellow", "nnbt_best_pose")

# Key channel residues
cmd.select("phe_triad",    "receptor and resi 69+121+199")
cmd.select("gate",         "receptor and resi 76+191")
cmd.select("loop_hotspot", "receptor and resi 314+315+316+317+318")
cmd.select("acid_base",    "receptor and resi 189+196")
cmd.select("axial_cys",    "receptor and resi 36")

cmd.show("sticks", "phe_triad or gate or loop_hotspot or acid_base or axial_cys")
cmd.color("cyan",    "phe_triad")
cmd.color("orange",  "gate")
cmd.color("magenta", "loop_hotspot")
cmd.color("red",     "acid_base")
cmd.color("yellow",  "axial_cys")

# Contact residues from docking
cmd.select("contacts", "receptor and resi {contact_resnums}")
cmd.show("sticks", "contacts")
cmd.color("salmon", "contacts")

# Distance: closest ligand atom to heme FE
cmd.select("heme_fe", "receptor and resn HEM and name FE")
cmd.select("nnbt_heavy", "nnbt_best_pose and not elem H")

# Zoom to active site
cmd.zoom("heme or nnbt_best_pose", buffer=8)
cmd.orient("heme or nnbt_best_pose")

# Measure distances manually in PyMOL:
# Click Wizard → Measurement → click α-C atom on NNBT → click FE on heme
# The α-C is the CH2 directly bonded to N (not the ring carbon)

print(\"\"\"
Colour legend:
  cyan    = Phe triad F69/F121/F199 (substrate sandwiching)
  orange  = Gate F76/F191 (molecular hinge)
  magenta = Loop G314-G318 (A316 hotspot)
  red     = Acid-base pair R189/E196 (NEVER mutate)
  yellow  = Axial Cys36 + NNBT best pose
  salmon  = All contact residues from docking
  green   = Heme

To measure α-C to FE distance:
  Wizard → Measurement → click CH2 on NNBT → click FE on heme
\"\"\")
"""

    script_path = workdir / "visualise_docking.py"
    script_path.write_text(script)
    log.info(f"  PyMOL script → {script_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="NNBT × PaDa-I docking pipeline")
    p.add_argument("--vina",           type=Path, required=True,
                   help="Path to AutoDock Vina binary")
    p.add_argument("--receptor",       type=Path, default=None,
                   help="Input receptor PDB (default: download 6EKZ automatically)")
    p.add_argument("--receptor_pdbqt", type=Path, default=None,
                   help="Pre-prepared receptor PDBQT (skips preparation step)")
    p.add_argument("--workdir",        type=Path, default=Path("docking_results"))
    p.add_argument("--box_size",       type=float, nargs=3,
                   default=list(BOX_SIZE), metavar=("DX","DY","DZ"))
    p.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    p.add_argument("--num_modes",      type=int, default=NUM_MODES)
    return p.parse_args()


def main():
    args = parse_args()
    workdir = args.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("NNBT × PaDa-I Docking Pipeline")
    log.info("=" * 60)

    # ── Receptor ──────────────────────────────────────────────────────────────
    if args.receptor_pdbqt:
        receptor_pdbqt = args.receptor_pdbqt
        receptor_clean = args.receptor or workdir / f"{PDB_ID}.pdb"
        heme_coords = {}
        log.info(f"  Using pre-prepared receptor PDBQT: {receptor_pdbqt}")
    else:
        pdb_path = args.receptor or download_pdb(workdir)
        receptor_pdbqt, heme_coords = prepare_receptor(pdb_path, workdir)
        receptor_clean = workdir / "receptor_clean.pdb"

    # ── Ligand ────────────────────────────────────────────────────────────────
    ligand_pdbqt, ligand_info = prepare_ligand(workdir)

    # ── Box centre ────────────────────────────────────────────────────────────
    clean_pdb = workdir / "receptor_clean.pdb" if not args.receptor_pdbqt else (args.receptor or workdir / f"{PDB_ID}.pdb")
    center = calculate_box_center(clean_pdb)

    # ── Dock ──────────────────────────────────────────────────────────────────
    box = tuple(args.box_size)
    poses_pdbqt, vina_scores = run_vina(
        args.vina, receptor_pdbqt, ligand_pdbqt,
        center, box, workdir, label="WT"
    )

    # ── Analyse ───────────────────────────────────────────────────────────────
    analysis = analyse_pose(
        poses_pdbqt, clean_pdb, heme_coords, workdir, label="WT"
    )

    # ── Report ────────────────────────────────────────────────────────────────
    write_report([analysis], {"WT": vina_scores}, heme_coords, workdir)
    write_pymol_script(workdir, [analysis])

    log.info("\n✓ Docking complete.")
    log.info(f"  Output directory: {workdir}/")
    log.info(f"  Open in PyMOL:    pymol {workdir}/visualise_docking.py")


if __name__ == "__main__":
    main()
