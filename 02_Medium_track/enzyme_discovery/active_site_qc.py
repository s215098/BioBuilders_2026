"""
active_site_qc.py
=================
Step 2c: After structure retrieval and classification, verify active site
integrity for each query structure.

Critical checks (cause pass_qc = False):
  1. Cys axial ligand — Cys SG within 3 Å of heme Fe (universal UPO feature)
  2. Heme present in experimental structures (predicted structures lack heme by
     design and are flagged informational-only)

Informational annotations (reported but do NOT affect pass/fail):
  3. Aromatic residues near Fe — counts Phe/Trp/Tyr within 8 Å of Fe.
     The Phe triad (Phe69/121/199 in AaeUPO) is well-characterised for
     PaDa-I/AaeUPO but conservation across other UPO clades is unknown.
     The raw count is reported for manual inspection rather than used as a
     pass/fail criterion.
  4. Channel residue count — total residues within 8 Å of Fe (fold sanity check)
  5. CA-RMSD to reference — fold-level similarity; high RMSD is expected for
     divergent clades and does not indicate a broken active site.

Requires: Biopython, numpy
Outputs:  results/active_site_qc.tsv
"""

import logging
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

from Bio.PDB import PDBParser, Superimposer, NeighborSearch
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Selection

logger = logging.getLogger(__name__)

PARSER = PDBParser(QUIET=True)


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class ActiveSiteQC:
    accession: str
    structure_source: str
    has_cys_axial: bool = False
    cys_axial_residue: Optional[str] = None    # e.g. "CYS_36_A"
    fe_coord: Optional[tuple] = None           # (x, y, z) of heme Fe
    n_residues_8A: int = 0                     # channel crowding proxy
    n_aromatics_near_fe: int = 0               # Phe/Trp/Tyr within 8 Å of Fe (informational)
    rmsd_to_ref: Optional[float] = None        # CA RMSD after superposition (informational)
    pass_qc: bool = False
    qc_flags: list = field(default_factory=list)  # CRITICAL flags cause pass_qc=False; INFO_ flags do not


# ---------------------------------------------------------------------------
# Core QC functions
# ---------------------------------------------------------------------------

def load_structure(path: Path, structure_id: str):
    """Load a PDB structure, return first model."""
    try:
        struct = PARSER.get_structure(structure_id, str(path))
        return struct[0]  # first model
    except Exception as e:
        logger.error(f"Cannot load structure {path}: {e}")
        return None


def find_heme_fe(model) -> Optional[np.ndarray]:
    """
    Locate heme iron (FE atom in HEM/HEC residue) in the structure.
    Returns numpy array [x, y, z] or None.
    """
    for chain in model:
        for res in chain:
            if res.resname in ("HEM", "HEC", "HEO", "FE2"):
                for atom in res:
                    if atom.element == "FE" or atom.name in ("FE", "FE2"):
                        return np.array(atom.coord)
    return None


def find_cys_axial(model, fe_coord: np.ndarray,
                   max_dist: float = 3.0) -> Optional[str]:
    """
    Find Cys residue whose S atom is within max_dist Å of heme Fe.
    This is the axial thiolate ligand diagnostic of UPO/CPO/P450.
    """
    for chain in model:
        for res in chain:
            if res.resname == "CYS":
                for atom in res:
                    if atom.name in ("SG", "S"):
                        dist = np.linalg.norm(atom.coord - fe_coord)
                        if dist <= max_dist:
                            return f"CYS_{res.id[1]}_{chain.id}"
    return None


def get_residues_near_fe(model, fe_coord: np.ndarray,
                         radius: float = 8.0) -> list:
    """
    Return list of residues with any heavy atom within radius Å of Fe.
    Used as a channel-crowding proxy.
    """
    atoms = list(Selection.unfold_entities(model, "A"))
    ns = NeighborSearch(atoms)
    nearby = ns.search(fe_coord, radius, "R")  # search for Residues
    return [r for r in nearby if r.resname not in ("HEM", "HEC", "HEO")]


AROMATIC_RESNAMES = {"PHE", "TRP", "TYR"}

def count_aromatics_near_fe(model, fe_coord: np.ndarray,
                             search_radius: float = 8.0) -> int:
    """
    Count aromatic residues (Phe/Trp/Tyr) within search_radius Å of Fe.
    Reported as an informational annotation — conservation of specific
    aromatic positions across UPO clades is not assumed.
    """
    nearby_res = get_residues_near_fe(model, fe_coord, radius=search_radius)
    return sum(1 for r in nearby_res if r.resname in AROMATIC_RESNAMES)


def compute_ca_rmsd(query_model, ref_model,
                    n_atoms: int = 50) -> Optional[float]:
    """
    Superimpose query onto reference using CA atoms, return RMSD.
    Uses the first n_atoms matched CA atoms.
    """
    try:
        query_cas = [a for a in Selection.unfold_entities(query_model, "A")
                     if a.name == "CA"][:n_atoms]
        ref_cas   = [a for a in Selection.unfold_entities(ref_model, "A")
                     if a.name == "CA"][:n_atoms]

        min_len = min(len(query_cas), len(ref_cas))
        if min_len < 20:
            return None

        sup = Superimposer()
        sup.set_atoms(ref_cas[:min_len], query_cas[:min_len])
        return round(sup.rms, 3)
    except Exception as e:
        logger.warning(f"RMSD calculation failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Per-structure QC
# ---------------------------------------------------------------------------

def run_qc_for_structure(accession: str, source: str,
                         structure_path: Optional[Path],
                         ref_model,
                         ref_phe_triad: Optional[list] = None) -> ActiveSiteQC:  # ref_phe_triad unused — kept for API compatibility
    """
    Full QC for one structure against a reference model.
    """
    qc = ActiveSiteQC(accession=accession, structure_source=source)

    if structure_path is None or not structure_path.exists():
        qc.qc_flags.append("NO_STRUCTURE")
        return qc

    model = load_structure(structure_path, accession)
    if model is None:
        qc.qc_flags.append("LOAD_FAILED")
        return qc

    # 1. Locate heme Fe
    fe_coord = find_heme_fe(model)
    if fe_coord is None:
        # Predicted structures (AlphaFold, Boltz) don't have heme — expected, not a failure
        if source in ("AlphaFold", "BOLTZ", "Boltz2"):
            qc.qc_flags.append("INFO_NO_HEME_PREDICTED_OK")
            logger.info(f"[{accession}] Predicted structure — heme must be grafted before docking")
            qc.rmsd_to_ref = compute_ca_rmsd(model, ref_model)
            qc.pass_qc = True
        else:
            qc.qc_flags.append("NO_HEME_EXPERIMENTAL")
            logger.warning(f"[{accession}] Experimental structure missing heme — check chain/residue")
            qc.rmsd_to_ref = compute_ca_rmsd(model, ref_model)
            qc.pass_qc = False
        return qc

    qc.fe_coord = tuple(float(x) for x in fe_coord)

    # 2. Cys axial ligand
    cys = find_cys_axial(model, fe_coord)
    if cys:
        qc.has_cys_axial = True
        qc.cys_axial_residue = cys
    else:
        qc.qc_flags.append("NO_CYS_AXIAL")
        logger.warning(f"[{accession}] No Cys axial ligand found within 3 Å of Fe")

    # 3. Channel residues (informational)
    channel_res = get_residues_near_fe(model, fe_coord, radius=8.0)
    qc.n_residues_8A = len(channel_res)
    if qc.n_residues_8A < 5:
        qc.qc_flags.append("INFO_SPARSE_CHANNEL")

    # 4. Aromatic residues near Fe (informational — conservation not assumed across clades)
    qc.n_aromatics_near_fe = count_aromatics_near_fe(model, fe_coord)
    if qc.n_aromatics_near_fe == 0:
        qc.qc_flags.append("INFO_NO_AROMATICS_NEAR_FE")
    logger.info(f"[{accession}] {qc.n_aromatics_near_fe} aromatic(s) (Phe/Trp/Tyr) within 8 Å of Fe")

    # 5. RMSD to reference (informational)
    qc.rmsd_to_ref = compute_ca_rmsd(model, ref_model)
    if qc.rmsd_to_ref and qc.rmsd_to_ref > 3.0:
        qc.qc_flags.append(f"INFO_HIGH_RMSD({qc.rmsd_to_ref})")

    # Pass/fail
    critical_flags = {"NO_CYS_AXIAL", "LOAD_FAILED", "NO_HEME_EXPERIMENTAL"}
    qc.pass_qc = not any(f.split("(")[0] in critical_flags for f in qc.qc_flags)

    logger.info(
        f"[{accession}] QC: Cys={qc.has_cys_axial} aromatics_near_fe={qc.n_aromatics_near_fe} "
        f"RMSD={qc.rmsd_to_ref} pass={qc.pass_qc} flags={qc.qc_flags}"
    )
    return qc


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_qc_results(qc_results: list[ActiveSiteQC], output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / "active_site_qc.tsv"

    with open(tsv_path, "w") as f:
        f.write(
            "accession\tsource\thas_cys_axial\tcys_residue\t"
            "fe_coord\tn_residues_8A\tn_aromatics_near_fe\t"
            "rmsd_to_ref\tpass_qc\tqc_flags\n"
        )
        for qc in qc_results:
            f.write(
                f"{qc.accession}\t{qc.structure_source}\t{qc.has_cys_axial}\t"
                f"{qc.cys_axial_residue or ''}\t"
                f"{qc.fe_coord or ''}\t{qc.n_residues_8A}\t"
                f"{qc.n_aromatics_near_fe}\t{qc.rmsd_to_ref or ''}\t"
                f"{qc.pass_qc}\t{';'.join(qc.qc_flags)}\n"
            )

    passed = sum(1 for q in qc_results if q.pass_qc)
    failed = len(qc_results) - passed
    logger.info(f"Active site QC → {tsv_path}  (pass:{passed} fail:{failed})")
    return tsv_path
