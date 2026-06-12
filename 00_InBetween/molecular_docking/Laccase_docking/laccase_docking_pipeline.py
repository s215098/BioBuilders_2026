#!/usr/bin/env python3
"""
Laccase Virtual Screening Pipeline
===================================
Docks two ligands (NNBT, Guaiacol) against ~79 laccase-like sequences using
AutoDock Vina at pH 7.  The T1 copper site is located by a single up-front
multiple sequence alignment (MSA) of the whole superfamily against the
reference: the alignment columns of the three copper-coordinating residues
(His443, Cys500, His505 in the reference) are mapped back to absolute residue
numbers in each target, whose sidechain atoms (ND1/NE2 for His, SG for Cys)
define the centroid of a 22 Å search box.  A structural Cα-superposition
fallback handles sequences too divergent to map through the MSA columns.

Each protein×ligand combination is docked 3× (seeds 1-3, exhaustiveness 24).
Mean and SD of best-pose affinity are reported alongside a geometry filter
that checks whether the ligand's redox-active atom lands within 8 Å of the
T1 centroid.  Receptors are protonated at pH 7 via pdb2pqr/PROPKA + meeko
before docking.

Dependencies (install once):
  pip install biopython rdkit requests tqdm pdb2pqr
  conda install -c conda-forge vina
  conda install -c conda-forge openbabel   # last-resort PDBQT fallback
  conda install -c bioconda mafft          # superfamily MSA (clustalo also OK)
  # If no MSA binary is present the pipeline falls back to a pure-Python
  # star alignment (Bio.Align.PairwiseAligner) anchored on the reference.

Directory layout expected / created:
  .
  ├── visible_sequences.csv          # input – Header, Sequence columns
  ├── reference.pdb                  # input – any laccase PDB with CU HETATM
  ├── itl_structures/                # input – PDB files for ItL01, ItL02, ItL03
  │   ├── ItL01.pdb
  │   ├── ItL02.pdb
  │   └── ItL03.pdb
  ├── structures/                    # auto-created – downloaded AlphaFold PDBs
  ├── trimmed/                       # auto-created – signal-peptide-trimmed PDBs
  ├── prepared/                      # auto-created – PDBQT receptor files
  ├── ligands/                       # auto-created – PDBQT ligand files
  └── results/                       # auto-created – docking poses + CSV summary

Usage:
  python laccase_docking_pipeline.py [options]

Options:
  --csv             Path to input CSV   (default: visible_sequences.csv)
  --reference       Path to reference PDB with T1 Cu (default: reference.pdb)
  --msa-tool        MSA engine: mafft | clustalo      (default: mafft)
  --itl-dir         Directory with ItL*.pdb files     (default: itl_structures)
  --box-size        Search box full-edge in Å          (default: 22)
  --exhaustiveness  Vina exhaustiveness per seed       (default: 24)
  --n-poses         Number of poses to save per run    (default: 9)
  --cpus            CPU threads for Vina               (default: all available)
  --dry-run         Print plan without running Vina
  --positive-control  Dock Guaiacol into reference.pdb before main loop
"""

import argparse
import csv
import json
import logging
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import List, Optional, Tuple

# ── optional heavy deps (imported lazily below) ───────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Ligand definitions
# ─────────────────────────────────────────────────────────────────────────────

LIGANDS = {
    "NNBT":     "CC1=CC=C(C=C1)N(CC(C)O)CC(C)O",
    "Guaiacol": "COC1=CC=CC=C1O",
}

# SMARTS for the redox-active atom in each ligand
# Guaiacol: phenolic OH; NNBT: tertiary amine nitrogen on the aromatic ring
LIGAND_REACTIVE_SMARTS = {
    "Guaiacol": "[OH1]c",
    "NNBT":     "[NX3;H0](c)",
}

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Dependency checks
# ─────────────────────────────────────────────────────────────────────────────

def check_dependencies() -> None:
    missing = []

    try:
        import Bio  # noqa: F401
    except ImportError:
        missing.append("biopython  →  pip install biopython")

    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import AllChem  # noqa: F401
    except ImportError:
        missing.append("rdkit      →  conda install -c conda-forge rdkit  OR  pip install rdkit")

    try:
        from vina import Vina  # noqa: F401
        log.info("Vina Python bindings found.")
    except ImportError:
        if not shutil.which("vina"):
            missing.append(
                "vina       →  conda install -c conda-forge vina  "
                "OR  download from https://github.com/ccsb-scripps/AutoDock-Vina/releases"
            )
        else:
            log.info("AutoDock Vina CLI found at: %s", shutil.which("vina"))

    # MSA engine is optional — a pure-Python fallback exists — so just inform.
    msa_bin = shutil.which("mafft") or shutil.which("clustalo")
    if msa_bin:
        log.info("MSA engine found at: %s", msa_bin)
    else:
        log.warning(
            "No MSA binary (mafft/clustalo) on PATH — T1 mapping will use the "
            "pure-Python star-alignment fallback. Install with: "
            "conda install -c bioconda mafft"
        )

    if missing:
        log.error("Missing required packages:\n  %s", "\n  ".join(missing))
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Structure acquisition
# ─────────────────────────────────────────────────────────────────────────────

ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uid}"

ITL_IDS = {"ItL01", "ItL02", "ItL03"}

def download_alphafold(uniprot_id: str, dest: Path) -> bool:
    """Download an AlphaFold model via the API; return True on success."""
    if dest.exists():
        log.debug("Already downloaded: %s", dest.name)
        return True
    try:
        log.info("Downloading %s …", uniprot_id)
        api_url = ALPHAFOLD_API_URL.format(uid=uniprot_id)
        req = urllib.request.Request(api_url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            payload = json.loads(r.read())
        pdb_url = payload[0]["pdbUrl"]
        with urllib.request.urlopen(pdb_url, timeout=60) as r, open(dest, "wb") as f:
            f.write(r.read())
        if dest.stat().st_size < 500:
            log.warning("  Download for %s looks truncated (%d bytes) — skipping.",
                        uniprot_id, dest.stat().st_size)
            dest.unlink()
            return False
        return True
    except Exception as exc:
        log.warning("  Failed to download %s: %s", uniprot_id, exc)
        if dest.exists():
            dest.unlink()
        return False


def collect_structures(
    df_rows: list,
    struct_dir: Path,
    itl_dir: Path,
) -> dict:
    """Returns dict {header: Path_to_pdb} for every successfully resolved entry."""
    struct_dir.mkdir(exist_ok=True)
    resolved = {}

    for row in df_rows:
        uid = row["Header"]

        if uid in ITL_IDS:
            candidate = itl_dir / f"{uid}.pdb"
            if candidate.exists():
                resolved[uid] = candidate
                log.info("ItL structure found: %s", candidate)
            else:
                log.warning("No PDB found for %s at %s – skipping.", uid, candidate)
            continue

        dest = struct_dir / f"{uid}.pdb"
        if download_alphafold(uid, dest):
            resolved[uid] = dest

    return resolved


# ─────────────────────────────────────────────────────────────────────────────
# 3.  T1 copper localisation
# ─────────────────────────────────────────────────────────────────────────────

# Sequence indices (0-based, in ca_and_seq output for reference.pdb chain A)
# and sidechain atoms to use for the T1 centroid.
# His443 → seq_idx 389, Cys500 → seq_idx 446, His505 → seq_idx 451
_REF_T1_RESIDUES: List[Tuple[int, Tuple[str, ...]]] = [
    (389, ("ND1", "NE2")),   # His443
    (446, ("SG",)),           # Cys500
    (451, ("ND1", "NE2")),   # His505
]

# Label / atom / expected-residue spec for the three T1 coordinating residues,
# in the same order as _REF_T1_RESIDUES.  Used by the MSA-mapping code paths.
_T1_SPEC: List[Tuple[str, Tuple[str, ...], str]] = [
    ("His1_num", ("ND1", "NE2"), "HIS"),  # His443
    ("Cys_num",  ("SG",),        "CYS"),  # Cys500
    ("His2_num", ("ND1", "NE2"), "HIS"),  # His505
]

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _ca_and_seq(struct) -> Tuple[str, list, list]:
    """Return (sequence_str, ca_atoms, residue_objects) for the first model."""
    seq_chars, cas, residues = [], [], []
    model = next(iter(struct))
    for chain in model:
        for residue in chain:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            seq_chars.append(_THREE_TO_ONE.get(residue.resname, "X"))
            cas.append(residue["CA"])
            residues.append(residue)
    return "".join(seq_chars), cas, residues


def _seq_from_pdb(pdb_path: Path) -> str:
    """One-letter sequence (Cα order, first model) of a PDB file."""
    from Bio.PDB import PDBParser  # type: ignore

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", str(pdb_path))
    seq, _, _ = _ca_and_seq(struct)
    return seq


def _build_aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices  # type: ignore
    aligner = PairwiseAligner()
    aligner.mode = "global"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        pass
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def get_cu_coords_from_pdb(pdb_path: Path) -> Optional[Tuple[float, float, float]]:
    """Return (x,y,z) of the *first* CU atom in HETATM records."""
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("HETATM") and line[12:16].strip() in ("CU", "CU1", "CU2"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                return (x, y, z)
    return None


# ── 3a.  Superfamily MSA + residue-number mapping ─────────────────────────────

REF_MSA_KEY = "__reference__"


def _parse_fasta(text: str) -> dict:
    """Parse aligned/unaligned FASTA text into an ordered {id: sequence} dict."""
    records: dict = {}
    header: Optional[str] = None
    chunks: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip().split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def _run_external_msa(sequences: dict, tool_bin: str, tool_name: str) -> dict:
    """Align `sequences` with an external mafft/clustalo binary; return {id: aln}."""
    with tempfile.TemporaryDirectory() as td:
        in_fa = Path(td) / "msa_in.fasta"
        out_fa = Path(td) / "msa_out.fasta"
        with open(in_fa, "w") as fh:
            for sid, seq in sequences.items():
                fh.write(f">{sid}\n{seq}\n")

        if tool_name == "mafft":
            cmd = [tool_bin, "--auto", "--quiet", str(in_fa)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip()[-300:])
            aligned_text = result.stdout
        elif tool_name == "clustalo":
            cmd = [tool_bin, "-i", str(in_fa), "-o", str(out_fa),
                   "--outfmt=fasta", "--force"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0 or not out_fa.exists():
                raise RuntimeError(result.stderr.strip()[-300:])
            aligned_text = out_fa.read_text()
        else:
            raise ValueError(f"Unknown MSA tool: {tool_name}")

        alignment = _parse_fasta(aligned_text)
        if set(alignment) != set(sequences):
            raise RuntimeError("MSA output is missing input sequences.")
        return alignment


def _gapped_pair(seq_a: str, seq_b: str, blocks) -> Tuple[str, str]:
    """Reconstruct the two gapped strings of a PairwiseAligner alignment."""
    a_blocks, b_blocks = blocks
    ra: List[str] = []
    rb: List[str] = []
    a_prev = b_prev = 0
    for (a_s, a_e), (b_s, b_e) in zip(a_blocks, b_blocks):
        # Unaligned residues before this matched block become mutual gaps.
        ra.append(seq_a[a_prev:a_s]); rb.append("-" * (a_s - a_prev))
        ra.append("-" * (b_s - b_prev)); rb.append(seq_b[b_prev:b_s])
        ra.append(seq_a[a_s:a_e]); rb.append(seq_b[b_s:b_e])
        a_prev, b_prev = a_e, b_e
    ra.append(seq_a[a_prev:]); rb.append("-" * (len(seq_a) - a_prev))
    ra.append("-" * (len(seq_b) - b_prev)); rb.append(seq_b[b_prev:])
    return "".join(ra), "".join(rb)


def _star_msa_fallback(sequences: dict, ref_id: str) -> dict:
    """
    Pure-Python star alignment anchored on the reference sequence.

    Each non-reference sequence is pairwise-aligned to the reference; the
    pairwise alignments are merged into one equal-length matrix by reserving,
    at every reference gap-slot, enough columns for the longest insertion seen.
    All residues of every sequence are preserved, so downstream residue-number
    counting stays correct.
    """
    aligner = _build_aligner()
    ref_seq = sequences[ref_id]
    n_ref = len(ref_seq)

    # Per target: matched[ref_pos] = aligned char (or '-'); ins[slot] = inserted run.
    per_target: dict = {}
    slot_max = [0] * (n_ref + 1)
    for sid, seq in sequences.items():
        if sid == ref_id:
            continue
        try:
            best = aligner.align(ref_seq, seq)[0]
            r_aln, t_aln = _gapped_pair(ref_seq, seq, best.aligned)
        except Exception as exc:
            log.warning("    Pairwise alignment failed for %s: %s — skipped.", sid, exc)
            continue

        matched = ["-"] * n_ref
        insertions = [""] * (n_ref + 1)
        ref_pos = -1
        for rc, tc in zip(r_aln, t_aln):
            if rc != "-":
                ref_pos += 1
                matched[ref_pos] = tc
            elif tc != "-":
                insertions[ref_pos + 1] += tc
        per_target[sid] = (matched, insertions)
        for i in range(n_ref + 1):
            slot_max[i] = max(slot_max[i], len(insertions[i]))

    def _build_row(matched: List[str], insertions: List[str]) -> str:
        out = [insertions[0].ljust(slot_max[0], "-")]
        for p in range(n_ref):
            out.append(matched[p])
            out.append(insertions[p + 1].ljust(slot_max[p + 1], "-"))
        return "".join(out)

    alignment = {ref_id: _build_row(list(ref_seq), [""] * (n_ref + 1))}
    for sid, (matched, insertions) in per_target.items():
        alignment[sid] = _build_row(matched, insertions)
    return alignment


def generate_superfamily_msa(
    sequences: dict,
    msa_tool: str = "mafft",
    ref_id: str = REF_MSA_KEY,
) -> dict:
    """
    Align the reference sequence with every target sequence into one MSA.

    Primary path uses an external `mafft`/`clustalo` binary via subprocess;
    if the binary is absent or fails, a pure-Python star alignment anchored on
    the reference is used instead (with a clear warning).  Returns an ordered
    {protein_id: gapped_aligned_sequence} dict with equal-length rows.
    """
    if len(sequences) < 2:
        return dict(sequences)

    tool_bin = shutil.which(msa_tool)
    if tool_bin:
        log.info("Running superfamily MSA with %s on %d sequences …",
                 msa_tool, len(sequences))
        try:
            alignment = _run_external_msa(sequences, tool_bin, msa_tool)
            n_cols = len(next(iter(alignment.values())))
            log.info("MSA complete via %s: %d sequences × %d columns.",
                     msa_tool, len(alignment), n_cols)
            return alignment
        except Exception as exc:
            log.warning("MSA tool '%s' failed (%s) — falling back to pure-Python "
                        "star alignment.", msa_tool, exc)
    else:
        log.warning("MSA tool '%s' not found on PATH — falling back to pure-Python "
                    "star alignment (Bio.Align.PairwiseAligner).", msa_tool)

    log.info("Building pure-Python star alignment on %d sequences …", len(sequences))
    alignment = _star_msa_fallback(sequences, ref_id=ref_id)
    n_cols = len(next(iter(alignment.values()))) if alignment else 0
    log.info("Star alignment complete: %d sequences × %d columns.",
             len(alignment), n_cols)
    return alignment


def map_t1_residues(alignment: dict, ref_id: str = REF_MSA_KEY) -> dict:
    """
    From an MSA, map the reference T1 columns to absolute residue numbers in
    every target sequence.

    Returns {protein_id: {"His1_num": X, "Cys_num": Y, "His2_num": Z}} where the
    numbers are 1-based positions in the target's *full* (ungapped) sequence —
    i.e. the residue numbering used by AlphaFold models.  A value is None if the
    target carries a gap at that alignment column.
    """
    if ref_id not in alignment:
        log.warning("Reference key %r absent from MSA — no T1 mapping produced.", ref_id)
        return {}

    ref_row = alignment[ref_id]

    # Reference ungapped-index → alignment column.
    col_of_ref_idx: dict = {}
    ungapped = -1
    for col, ch in enumerate(ref_row):
        if ch != "-":
            ungapped += 1
            col_of_ref_idx[ungapped] = col

    labels = [label for label, _, _ in _T1_SPEC]
    t1_cols: List[int] = []
    for (ref_idx, _), label in zip(_REF_T1_RESIDUES, labels):
        if ref_idx not in col_of_ref_idx:
            log.warning("Reference T1 residue at seq-index %d is outside the MSA "
                        "(%s) — mapping aborted.", ref_idx, label)
            return {}
        t1_cols.append(col_of_ref_idx[ref_idx])

    log.info("T1 reference columns located: %s → MSA columns %s",
             labels, t1_cols)

    # Expected one-letter residue at each T1 column (HIS→H, CYS→C, …).
    expected = [_THREE_TO_ONE.get(expect, "X") for _, _, expect in _T1_SPEC]

    mapping: dict = {}
    flagged: List[Tuple[str, List[str]]] = []
    for pid, row in alignment.items():
        if pid == ref_id:
            continue
        resnums: dict = {}
        bad: List[str] = []
        for label, col, exp in zip(labels, t1_cols, expected):
            if col >= len(row) or row[col] == "-":
                resnums[label] = None
                bad.append(f"{label}=gap")
            else:
                resnums[label] = sum(1 for ch in row[: col + 1] if ch != "-")
                if row[col] != exp:
                    bad.append(f"{label}={row[col]}≠{exp}")
        mapping[pid] = resnums
        if bad:
            flagged.append((pid, bad))

    n_ok = len(mapping) - len(flagged)
    log.info("T1 columns conserved (His/Cys/His) in %d / %d targets.",
             n_ok, len(mapping))
    if flagged:
        log.warning(
            "%d target(s) do NOT present His/Cys/His at the T1 columns and will "
            "use the structural superposition fallback for localisation:",
            len(flagged),
        )
        for pid, bad in sorted(flagged):
            log.warning("    %-14s %s", pid, ", ".join(bad))
    return mapping


# ── 3b.  Per-structure T1 centroid ────────────────────────────────────────────

def _t1_centroid_from_mapping(
    mobile_pdb: Path,
    t1_mapping: dict,
):
    """
    Centroid of the T1 sidechain atoms using pre-computed residue numbers.

    Returns the (x,y,z) centroid, or None if no expected residue could be
    validated (caller then drops to the superposition fallback).
    """
    from Bio.PDB import PDBParser  # type: ignore
    import numpy as np  # type: ignore

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("mob", str(mobile_pdb))
    model = next(iter(struct))

    by_num: dict = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            by_num.setdefault(res.id[1], res)   # first chain wins on clashes

    coords: List["np.ndarray"] = []
    n_ok = 0
    for label, atom_names, expect in _T1_SPEC:
        resnum = t1_mapping.get(label)
        if resnum is None:
            log.warning("    %s: no MSA column mapping for %s.", mobile_pdb.name, label)
            continue
        res = by_num.get(resnum)
        if res is None:
            log.warning("    %s: residue %d (%s) not present in PDB.",
                        mobile_pdb.name, resnum, label)
            continue
        if res.resname != expect:
            log.warning("    %s: residue %d is %s, expected %s (%s) — skipped.",
                        mobile_pdb.name, resnum, res.resname, expect, label)
            continue
        n_ok += 1
        for aname in atom_names:
            if aname in res:
                coords.append(res[aname].get_vector().get_array())
            else:
                log.warning("    %s: atom %s missing from residue %d (%s).",
                            mobile_pdb.name, aname, resnum, label)

    if not coords:
        return None

    if n_ok < len(_T1_SPEC):
        log.warning("    %s: only %d/%d T1 residues validated via MSA — partial centroid.",
                    mobile_pdb.name, n_ok, len(_T1_SPEC))

    centroid = np.mean(coords, axis=0)
    log.info("    MSA-mapped T1 centroid (%d/%d residues): (%.2f, %.2f, %.2f)",
             n_ok, len(_T1_SPEC), *centroid)
    return tuple(float(v) for v in centroid)


def locate_t1_site(
    mobile_pdb: Path,
    reference_pdb: Path,
    t1_mapping: Optional[dict] = None,
) -> Optional[Tuple[float, float, float]]:
    """
    Locate the T1 Cu site in mobile_pdb.

    Primary path uses the pre-computed MSA residue-number mapping
    (`t1_mapping = {"His1_num": X, "Cys_num": Y, "His2_num": Z}`): the sidechain
    atoms of those residues (ND1/NE2 for His, SG for Cys) are read straight from
    the PDB and averaged into a centroid.

    If no mapping is supplied, or the mapped residues cannot be validated (a
    sequence too divergent for the MSA columns to be meaningful), it falls back
    to the structural Cα-superposition route in `_locate_t1_via_superposition`.
    """
    if t1_mapping:
        centroid = _t1_centroid_from_mapping(mobile_pdb, t1_mapping)
        if centroid is not None:
            return centroid
        log.warning("  MSA mapping unusable for %s — using superposition fallback.",
                    mobile_pdb.name)
    else:
        log.info("  No MSA mapping for %s — using superposition fallback.",
                 mobile_pdb.name)

    return _locate_t1_via_superposition(mobile_pdb, reference_pdb)


def _locate_t1_via_superposition(
    mobile_pdb: Path,
    reference_pdb: Path,
) -> Optional[Tuple[float, float, float]]:
    """
    Structural fallback: map reference coordinating residues (His443, Cys500,
    His505) via on-the-fly sequence-guided Cα alignment + Superimposer, and
    return the centroid of their sidechain atoms.  Falls back further to
    superposition-based Cu-coordinate transfer if no residues map.
    """
    from Bio.PDB import PDBParser, Superimposer  # type: ignore
    import numpy as np  # type: ignore

    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", str(reference_pdb))
    mob_struct = parser.get_structure("mob", str(mobile_pdb))

    ref_seq, ref_cas, ref_residues = _ca_and_seq(ref_struct)
    mob_seq, mob_cas, mob_residues = _ca_and_seq(mob_struct)

    if len(ref_seq) < 10 or len(mob_seq) < 10:
        log.warning(
            "  Too few Cα atoms (%s vs ref). Falling back to raw ref Cu coords.",
            mobile_pdb.name,
        )
        return get_cu_coords_from_pdb(reference_pdb)

    aligner = _build_aligner()
    best = aligner.align(ref_seq, mob_seq)[0]

    ref_to_mob: dict = {}
    ref_paired, mob_paired = [], []
    for (r_s, r_e), (m_s, m_e) in zip(best.aligned[0], best.aligned[1]):
        for r_i, m_i in zip(range(r_s, r_e), range(m_s, m_e)):
            ref_to_mob[r_i] = m_i
            ref_paired.append(ref_cas[r_i])
            mob_paired.append(mob_cas[m_i])

    n_pairs = len(ref_paired)

    # Always compute superposition (needed for RMSD logging and fallback)
    sup = None
    rot = tran = None
    if n_pairs >= 10:
        sup = Superimposer()
        sup.set_atoms(ref_paired, mob_paired)
        rot, tran = sup.rotran

    # Map T1 coordinating residues → target sidechain atoms
    centroid_coords: List[np.ndarray] = []
    n_attempted = sum(len(atoms) for _, atoms in _REF_T1_RESIDUES)

    for ref_idx, atom_names in _REF_T1_RESIDUES:
        if ref_idx not in ref_to_mob:
            log.warning(
                "  T1 ref residue at seq-pos %d has no aligned partner in %s — skipped.",
                ref_idx, mobile_pdb.name,
            )
            continue
        mob_idx = ref_to_mob[ref_idx]
        mob_res = mob_residues[mob_idx]
        for aname in atom_names:
            if aname in mob_res:
                centroid_coords.append(mob_res[aname].get_vector().get_array())
            else:
                log.warning(
                    "  Atom %s missing from %s residue at pos %d — skipped.",
                    aname, mobile_pdb.name, mob_idx,
                )

    n_mapped = len(centroid_coords)

    if n_mapped == 0:
        # Full fallback: superposition-based Cu coordinate transfer
        log.warning(
            "  No T1 residue atoms mapped for %s — falling back to superposition Cu transfer.",
            mobile_pdb.name,
        )
        if rot is None:
            return get_cu_coords_from_pdb(reference_pdb)
        ref_cu = get_cu_coords_from_pdb(reference_pdb)
        if ref_cu is None:
            return None
        mob_cu = rot.T @ (np.array(ref_cu) - tran)
        rmsd_str = f"{sup.rms:.2f} Å" if sup else "?"
        log.info(
            "  Fallback: RMSD=%s (%d pairs)  →  T1 Cu=(%.1f, %.1f, %.1f)",
            rmsd_str, n_pairs, *mob_cu,
        )
        return tuple(float(v) for v in mob_cu)

    if n_mapped < n_attempted:
        log.warning(
            "  Only %d/%d T1 atoms mapped for %s — using partial centroid.",
            n_mapped, n_attempted, mobile_pdb.name,
        )

    centroid = np.mean(centroid_coords, axis=0)
    rmsd_str = f"{sup.rms:.2f} Å" if sup else "?"
    log.info(
        "  Alignment RMSD=%s (%d pairs), %d/%d T1 atoms  →  centroid=(%.2f, %.2f, %.2f)",
        rmsd_str, n_pairs, n_mapped, n_attempted, *centroid,
    )
    return tuple(float(v) for v in centroid)


# ─────────────────────────────────────────────────────────────────────────────
# 3c.  Signal-peptide / disordered N-terminal trim  (AlphaFold only)
# ─────────────────────────────────────────────────────────────────────────────

PLDDT_THRESHOLD = 70.0   # Cα B-factor below this is considered low-confidence


def trim_n_terminal(
    pdb_path: Path,
    reference_pdb: Path,
    trimmed_dir: Path,
    is_alphafold: bool = True,
) -> Tuple[Path, int]:
    """
    Strip N-terminal residues that are BOTH outside the reference-aligned core
    (i.e., N-terminal to the first aligned residue) AND have pLDDT < 70.

    Only applied to AlphaFold structures; ItL structures are returned unchanged.
    Returns (output_path, n_trimmed).
    """
    trimmed_dir.mkdir(exist_ok=True)

    if not is_alphafold:
        return pdb_path, 0

    trimmed_path = trimmed_dir / pdb_path.name

    try:
        from Bio.PDB import PDBParser, PDBIO, Select  # type: ignore

        parser = PDBParser(QUIET=True)
        ref_struct = parser.get_structure("ref", str(reference_pdb))
        mob_struct = parser.get_structure("mob", str(pdb_path))

        ref_seq, ref_cas, _ = _ca_and_seq(ref_struct)
        mob_seq, mob_cas, mob_residues = _ca_and_seq(mob_struct)

        if len(ref_seq) < 10 or len(mob_seq) < 10:
            return pdb_path, 0

        aligner = _build_aligner()
        best = aligner.align(ref_seq, mob_seq)[0]

        # First aligned position in mobile
        if best.aligned[1].size == 0:
            return pdb_path, 0
        first_mob_aligned = int(best.aligned[1][0][0])

        # Residues N-terminal to first_mob_aligned with low pLDDT
        to_trim_ids = set()
        for i in range(first_mob_aligned):
            res = mob_residues[i]
            ca = res["CA"]
            if ca.get_bfactor() < PLDDT_THRESHOLD:
                to_trim_ids.add((res.get_parent().id, res.id))

        if not to_trim_ids:
            return pdb_path, 0

        class _TrimSelect(Select):
            def accept_residue(self, res):
                key = (res.get_parent().id, res.id)
                return 0 if key in to_trim_ids else 1

        io = PDBIO()
        io.set_structure(mob_struct)
        io.save(str(trimmed_path), select=_TrimSelect())

        n = len(to_trim_ids)
        log.info("  Trimmed %d low-pLDDT N-terminal residues from %s → %s",
                 n, pdb_path.name, trimmed_path.name)
        return trimmed_path, n

    except Exception as exc:
        log.warning("  N-terminal trim failed for %s: %s — using original.", pdb_path.name, exc)
        return pdb_path, 0


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Receptor preparation  (PDB → PDBQT)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_receptor(pdb_path: Path, out_pdbqt: Path) -> bool:
    """
    Convert a protein PDB to PDBQT for Vina.  Preferred path:
      a) pdb2pqr (PROPKA, pH 7) → protonated PDB → meeko  [best quality]
      b) meeko directly on input PDB
      c) prepare_receptor4.py (MGLTools)
      d) obabel  [last resort — warn clearly]
    Returns True on success.
    """
    if out_pdbqt.exists():
        return True

    # Strip HETATM lines for pdb2pqr (metals confuse it)
    def _write_protein_only(src: Path, dst: Path) -> None:
        with open(src) as fi, open(dst, "w") as fo:
            for line in fi:
                if not line.startswith("HETATM"):
                    fo.write(line)

    # ── a) pdb2pqr + meeko ────────────────────────────────────────────────────
    pdb2pqr_bin = shutil.which("pdb2pqr") or shutil.which("pdb2pqr3") or shutil.which("pdb2pqr30")
    if pdb2pqr_bin:
        with tempfile.TemporaryDirectory() as td:
            clean_pdb  = Path(td) / "clean.pdb"
            prot_pdb   = Path(td) / "protonated.pdb"
            tmp_pqr    = Path(td) / "out.pqr"
            _write_protein_only(pdb_path, clean_pdb)
            cmd = [
                pdb2pqr_bin,
                "--ff", "PARSE",
                "--with-ph", "7.0",
                "--pdb-output", str(prot_pdb),
                str(clean_pdb), str(tmp_pqr),
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0 and prot_pdb.exists():
                # Convert protonated PDB → PDBQT with meeko
                try:
                    from meeko import PDBQTWriterLegacy  # type: ignore
                    from Bio.PDB import PDBParser  # type: ignore

                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure("rec", str(prot_pdb))
                    for model in structure:
                        for chain in model:
                            for res in list(chain):
                                if res.id[0] != " ":
                                    chain.detach_child(res.id)
                    writer = PDBQTWriterLegacy()
                    pdbqt_str = writer.write_string(structure)[0]
                    out_pdbqt.write_text(pdbqt_str)
                    log.info("  Receptor prepared via pdb2pqr+meeko (pH 7): %s", out_pdbqt.name)
                    return True
                except Exception as e:
                    log.debug("  meeko step after pdb2pqr failed: %s", e)
            elif pdb2pqr_bin:
                log.debug("  pdb2pqr failed: %s", result.stderr.strip()[:200])

    # ── b) meeko directly ─────────────────────────────────────────────────────
    try:
        from meeko import PDBQTWriterLegacy  # type: ignore
        from Bio.PDB import PDBParser  # type: ignore

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("rec", str(pdb_path))
        for model in structure:
            for chain in model:
                for res in list(chain):
                    if res.id[0] != " ":
                        chain.detach_child(res.id)

        writer = PDBQTWriterLegacy()
        pdbqt_str = writer.write_string(structure)[0]
        out_pdbqt.write_text(pdbqt_str)
        log.info("  Receptor prepared via meeko: %s", out_pdbqt.name)
        return True
    except Exception:
        pass

    # ── c) MGLTools prepare_receptor4.py ──────────────────────────────────────
    prep4 = shutil.which("prepare_receptor4.py") or shutil.which("prepare_receptor4")
    if prep4:
        cmd = [prep4, "-r", str(pdb_path), "-o", str(out_pdbqt),
               "-A", "hydrogens", "-U", "nphs_lps_waters_nonstdres"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and out_pdbqt.exists():
            log.info("  Receptor prepared via MGLTools: %s", out_pdbqt.name)
            return True
        log.warning("  MGLTools failed: %s", result.stderr.strip())

    # ── d) obabel last resort ─────────────────────────────────────────────────
    if shutil.which("obabel"):
        cmd = ["obabel", str(pdb_path), "-O", str(out_pdbqt), "-xr"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and out_pdbqt.exists():
            log.warning(
                "  Receptor prepared via obabel (low quality — no pH-7 protonation): %s",
                out_pdbqt.name,
            )
            return True
        log.warning("  obabel failed: %s", result.stderr.strip())

    log.error(
        "  Could not prepare receptor %s. "
        "Install meeko (pip install meeko) or pdb2pqr (pip install pdb2pqr).",
        pdb_path.name,
    )
    return False


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Ligand preparation  (SMILES → PDBQT + SDF)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ligand(name: str, smiles: str, out_pdbqt: Path) -> bool:
    """
    SMILES → 3-D conformer → PDBQT (via meeko or obabel).
    Also saves an SDF alongside the PDBQT for reactive-atom identification.
    """
    sdf_path = out_pdbqt.with_suffix(".sdf")
    if out_pdbqt.exists() and sdf_path.exists():
        return True

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log.error("  Invalid SMILES for %s: %s", name, smiles)
            return False

        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            log.error("  Conformer generation failed for %s.", name)
            return False
        AllChem.MMFFOptimizeMolecule(mol)

        # Save heavy-atom SDF for later reactive-atom identification
        mol_no_h = Chem.RemoveHs(mol)
        try:
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(mol_no_h)
            writer.close()
        except Exception as e:
            log.debug("  Could not save SDF for %s: %s", name, e)

        # ── meeko path ───────────────────────────────────────────────────────
        try:
            from meeko import MoleculePreparation  # type: ignore

            prep = MoleculePreparation()
            prep.prepare(mol)
            pdbqt_str = prep.write_pdbqt_string()
            out_pdbqt.write_text(pdbqt_str)
            log.info("  Ligand %s prepared via meeko.", name)
            return True
        except Exception as e:
            log.debug("  meeko ligand prep failed (%s), trying obabel.", e)

        # ── obabel fallback ──────────────────────────────────────────────────
        if shutil.which("obabel"):
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
                tmp_sdf = Path(tmp.name)
            writer2 = Chem.SDWriter(str(tmp_sdf))
            writer2.write(mol)
            writer2.close()
            cmd = ["obabel", str(tmp_sdf), "-O", str(out_pdbqt)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            tmp_sdf.unlink(missing_ok=True)
            if result.returncode == 0 and out_pdbqt.exists():
                log.info("  Ligand %s prepared via obabel.", name)
                return True
            log.warning("  obabel ligand failed: %s", result.stderr.strip())

        log.error(
            "  Cannot write PDBQT for %s. Install meeko or openbabel.", name,
        )
        return False

    except ImportError:
        log.error("  RDKit not installed. pip install rdkit")
        return False


# ─────────────────────────────────────────────────────────────────────────────
# 5b.  Reactive-atom helpers for T1 geometry filter
# ─────────────────────────────────────────────────────────────────────────────

def _read_pdbqt_atoms(
    pdbqt_path: Path,
    pose: int = 0,
) -> List[Tuple[str, float, float, float]]:
    """
    Return [(atom_name, x, y, z), ...] for heavy atoms in the given pose (0-based).
    Handles both multi-MODEL PDBQT (Vina output) and single-pose PDBQT (ligand input).
    """
    atoms: List[Tuple[str, float, float, float]] = []
    current_pose = -1
    cur: List[Tuple[str, float, float, float]] = []
    found = False

    try:
        with open(pdbqt_path) as fh:
            for line in fh:
                if line.startswith("MODEL"):
                    current_pose += 1
                    cur = []
                elif line.startswith("ENDMDL"):
                    if current_pose == pose:
                        atoms = cur
                        found = True
                        break
                elif (line.startswith("ATOM") or line.startswith("HETATM")) and current_pose == pose:
                    aname = line[12:16].strip()
                    if aname.startswith("H"):
                        continue
                    try:
                        cur.append((aname, float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except ValueError:
                        pass

        # Single-pose PDBQT (no MODEL/ENDMDL)
        if not found and current_pose == -1:
            with open(pdbqt_path) as fh:
                for line in fh:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        aname = line[12:16].strip()
                        if aname.startswith("H"):
                            continue
                        try:
                            atoms.append((aname, float(line[30:38]),
                                          float(line[38:46]), float(line[46:54])))
                        except ValueError:
                            pass
    except Exception:
        pass
    return atoms


def find_reactive_atom_name(ligand_pdbqt: Path, smarts: str) -> Optional[str]:
    """
    Return the PDBQT atom name of the reactive atom identified by SMARTS.
    Uses the companion SDF (saved by prepare_ligand) to get 3-D coordinates,
    then finds the closest heavy atom in the ligand PDBQT.
    """
    try:
        from rdkit import Chem
        import numpy as np

        sdf_path = ligand_pdbqt.with_suffix(".sdf")
        if not sdf_path.exists():
            return None

        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True)
        if mol is None:
            return None
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            return None
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            return None

        conf = mol.GetConformer()
        pos = conf.GetAtomPosition(matches[0][0])
        reactive_coord = np.array([pos.x, pos.y, pos.z])

        pdbqt_atoms = _read_pdbqt_atoms(ligand_pdbqt, pose=0)
        if not pdbqt_atoms:
            return None

        best_name, best_dist = None, float("inf")
        for aname, ax, ay, az in pdbqt_atoms:
            d = float(np.linalg.norm(np.array([ax, ay, az]) - reactive_coord))
            if d < best_dist:
                best_dist = d
                best_name = aname

        return best_name

    except Exception as e:
        log.debug("  find_reactive_atom_name failed: %s", e)
        return None


def measure_reactive_dist(
    pose_pdbqt: Path,
    reactive_atom_name: str,
    center: Tuple[float, float, float],
    pose_idx: int = 0,
) -> Optional[float]:
    """Return distance (Å) from the reactive atom in pose_idx to center, or None."""
    try:
        import numpy as np
        cx, cy, cz = center
        for aname, x, y, z in _read_pdbqt_atoms(pose_pdbqt, pose=pose_idx):
            if aname == reactive_atom_name:
                return float(math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2))
        return None
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Docking
# ─────────────────────────────────────────────────────────────────────────────

_REPLICATE_SEEDS = [1, 2, 3]
_REPLICATE_EXHAUSTIVENESS = 24


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: Tuple[float, float, float],
    box_size: float,
    out_pdbqt: Path,
    out_log: Path,
    exhaustiveness: int,
    n_poses: int,
    cpus: int,
    seed: Optional[int] = None,
) -> Optional[List[Tuple[int, float]]]:
    """
    Run AutoDock Vina (Python API or CLI).
    box_size is the full edge length in Å.
    Returns list of (pose_index, affinity_kcal_mol) or None on failure.
    """
    cx, cy, cz = center
    sz = box_size  # full edge — no multiplication

    # ── Python API (preferred) ──────────────────────────────────────────────
    try:
        from vina import Vina as VinaAPI  # type: ignore

        v = VinaAPI(sf_name="vina", cpu=cpus, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.compute_vina_maps(center=[cx, cy, cz], box_size=[sz, sz, sz])
        try:
            v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed)
        except TypeError:
            v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(str(out_pdbqt), n_poses=n_poses, overwrite=True)

        energies = v.energies(n_poses=n_poses)
        results = [(i + 1, float(e[0])) for i, e in enumerate(energies)]

        with open(out_log, "w") as fh:
            fh.write("pose,affinity_kcal_mol\n")
            for pose, aff in results:
                fh.write(f"{pose},{aff:.3f}\n")
        return results

    except ImportError:
        pass

    # ── CLI fallback ────────────────────────────────────────────────────────
    vina_bin = shutil.which("vina")
    if not vina_bin:
        log.error("AutoDock Vina not found.  Install: conda install -c conda-forge vina")
        return None

    cmd = [
        vina_bin,
        "--receptor",    str(receptor_pdbqt),
        "--ligand",      str(ligand_pdbqt),
        "--center_x",    f"{cx:.3f}",
        "--center_y",    f"{cy:.3f}",
        "--center_z",    f"{cz:.3f}",
        "--size_x",      f"{sz:.1f}",
        "--size_y",      f"{sz:.1f}",
        "--size_z",      f"{sz:.1f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",   str(n_poses),
        "--out",         str(out_pdbqt),
        "--log",         str(out_log),
        "--cpu",         str(cpus),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("  Vina CLI error:\n%s", result.stderr[-800:])
        return None

    affinities = []
    for line in result.stdout.splitlines() + list(open(out_log)):
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                affinities.append((int(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return affinities if affinities else []


def run_vina_replicates(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: Tuple[float, float, float],
    box_size: float,
    out_pdbqt: Path,
    out_log: Path,
    n_poses: int,
    cpus: int,
) -> Tuple[Optional[List[Tuple[int, float]]], float, float]:
    """
    Run 3 replicate docking runs (seeds 1, 2, 3 at exhaustiveness 24).
    Poses from seed 1 are written to out_pdbqt.
    Returns (seed1_results, best_affinity_mean, best_affinity_sd).
    """
    best_per_seed: List[float] = []
    seed1_results: Optional[List[Tuple[int, float]]] = None

    for seed in _REPLICATE_SEEDS:
        seed_pdbqt = out_pdbqt if seed == 1 else out_pdbqt.with_suffix(f".seed{seed}.pdbqt")
        seed_log   = out_log   if seed == 1 else out_log.with_suffix(f".seed{seed}.log")

        results = run_vina(
            receptor_pdbqt, ligand_pdbqt, center, box_size,
            seed_pdbqt, seed_log,
            exhaustiveness=_REPLICATE_EXHAUSTIVENESS,
            n_poses=n_poses,
            cpus=cpus,
            seed=seed,
        )
        if results:
            best_per_seed.append(min(a for _, a in results))
            if seed == 1:
                seed1_results = results

    if not best_per_seed:
        return None, float("nan"), float("nan")

    mean_aff = statistics.mean(best_per_seed)
    sd_aff   = statistics.stdev(best_per_seed) if len(best_per_seed) > 1 else 0.0
    return seed1_results, mean_aff, sd_aff


# ─────────────────────────────────────────────────────────────────────────────
# 6b.  Positive control
# ─────────────────────────────────────────────────────────────────────────────

def run_positive_control(
    ref_pdb: Path,
    guaiacol_pdbqt: Path,
    prep_dir: Path,
    res_dir: Path,
    box_size: float,
    exhaustiveness: int,
    n_poses: int,
    cpus: int,
) -> bool:
    """
    Dock Guaiacol into reference.pdb at the reference T1 Cu coordinates.
    PASS if the phenolic O of the top pose lands within 8 Å of the T1 Cu.
    Returns True on PASS.
    """
    log.info("═" * 60)
    log.info("POSITIVE CONTROL: Guaiacol → %s", ref_pdb.name)

    ref_cu = get_cu_coords_from_pdb(ref_pdb)
    if ref_cu is None:
        log.error("  Positive control: no CU found in reference — SKIPPED.")
        return False

    ref_pdbqt = prep_dir / f"_reference_ctrl.pdbqt"
    if not prepare_receptor(ref_pdb, ref_pdbqt):
        log.error("  Positive control: receptor prep failed — SKIPPED.")
        return False

    out_pdbqt = res_dir / "_positive_control_Guaiacol_poses.pdbqt"
    out_log   = res_dir / "_positive_control_Guaiacol.log"

    results = run_vina(
        ref_pdbqt, guaiacol_pdbqt, ref_cu, box_size,
        out_pdbqt, out_log,
        exhaustiveness=exhaustiveness,
        n_poses=n_poses,
        cpus=cpus,
        seed=1,
    )

    if not results or not out_pdbqt.exists():
        log.warning("  Positive control: docking failed — SKIPPED.")
        return False

    reactive_name = find_reactive_atom_name(guaiacol_pdbqt, LIGAND_REACTIVE_SMARTS["Guaiacol"])
    dist = None
    if reactive_name:
        dist = measure_reactive_dist(out_pdbqt, reactive_name, ref_cu, pose_idx=0)

    best_aff = min(a for _, a in results)
    if dist is not None:
        passed = dist <= 8.0
        status = "PASS" if passed else "FAIL"
        log.info(
            "  Positive control: best=%.2f kcal/mol, phenolic-O dist=%.1f Å → %s",
            best_aff, dist, status,
        )
        return passed
    else:
        log.warning(
            "  Positive control: best=%.2f kcal/mol; could not measure phenolic-O distance.",
            best_aff,
        )
        return False


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Summary CSV writer
# ─────────────────────────────────────────────────────────────────────────────

def write_summary(records: list, out_csv: Path) -> None:
    fieldnames = [
        "protein_id", "ligand",
        "best_affinity_mean", "best_affinity_sd", "best_affinity_kcal_mol",
        "all_affinities_kcal_mol",
        "t1_cu_x", "t1_cu_y", "t1_cu_z",
        "reactive_atom_dist_to_T1", "productive_pose", "ranking_score",
        "poses_pdbqt",
    ]
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    log.info("Summary written → %s", out_csv)


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Main orchestration
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv",            default="visible_sequences.csv")
    p.add_argument("--reference",      default="reference.pdb",
                   help="PDB file of a reference laccase with a CU HETATM record")
    p.add_argument("--msa-tool",       default="mafft", dest="msa_tool",
                   choices=["mafft", "clustalo"],
                   help="External MSA engine (default: mafft); pure-Python "
                        "fallback used if the binary is unavailable")
    p.add_argument("--itl-dir",        default="itl_structures", dest="itl_dir")
    p.add_argument("--box-size",       type=float, default=22.0,
                   help="Search box full-edge in Å (default: 22)")
    p.add_argument("--exhaustiveness", type=int,   default=_REPLICATE_EXHAUSTIVENESS,
                   help="Vina exhaustiveness per replicate seed (default: 24)")
    p.add_argument("--n-poses",        type=int,   default=9)
    p.add_argument("--cpus",           type=int,   default=os.cpu_count() or 4)
    p.add_argument("--dry-run",        action="store_true")
    p.add_argument("--positive-control", action="store_true", dest="positive_control",
                   help="Dock Guaiacol into reference.pdb before main loop as a sanity check")
    p.add_argument("--struct-dir",     default="structures")
    p.add_argument("--prepared-dir",   default="prepared")
    p.add_argument("--trimmed-dir",    default="trimmed")
    p.add_argument("--ligand-dir",     default="ligands")
    p.add_argument("--results-dir",    default="results")
    return p.parse_args()


def main():
    args = parse_args()

    check_dependencies()

    # ── Paths ──────────────────────────────────────────────────────────────
    csv_path    = Path(args.csv)
    ref_pdb     = Path(args.reference)
    itl_dir     = Path(args.itl_dir)
    struct_dir  = Path(args.struct_dir)
    trimmed_dir = Path(args.trimmed_dir)
    prep_dir    = Path(args.prepared_dir)
    lig_dir     = Path(args.ligand_dir)
    res_dir     = Path(args.results_dir)

    for d in (struct_dir, trimmed_dir, prep_dir, lig_dir, res_dir):
        d.mkdir(exist_ok=True)

    # ── Validate reference PDB ─────────────────────────────────────────────
    if not ref_pdb.exists():
        log.error("Reference PDB not found: %s", ref_pdb)
        sys.exit(1)

    ref_cu = get_cu_coords_from_pdb(ref_pdb)
    if ref_cu is None:
        log.error("No CU HETATM atom found in %s.", ref_pdb)
        sys.exit(1)
    log.info("Reference T1 Cu coordinates: (%.2f, %.2f, %.2f)", *ref_cu)

    # ── Read sequences CSV ─────────────────────────────────────────────────
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    log.info("Loaded %d sequences from %s", len(rows), csv_path)

    seq_by_header = {row["Header"]: row["Sequence"].strip() for row in rows}

    # ── Download / locate structures ────────────────────────────────────────
    structures = collect_structures(rows, struct_dir, itl_dir)
    log.info("%d / %d structures resolved.", len(structures), len(rows))

    # ── Pre-docking MSA phase: map T1 residues across the whole superfamily ──
    log.info("═" * 60)
    log.info("MSA PHASE: locating T1 coordinating residues across the superfamily")
    ref_seq = _seq_from_pdb(ref_pdb)
    msa_input = {REF_MSA_KEY: ref_seq}
    for uid in structures:
        seq = seq_by_header.get(uid)
        if seq:
            msa_input[uid] = seq
        else:
            log.warning("  No CSV sequence for %s — it will use the superposition "
                        "fallback for T1 localisation.", uid)

    alignment = generate_superfamily_msa(msa_input, msa_tool=args.msa_tool)
    t1_mapping = map_t1_residues(alignment, ref_id=REF_MSA_KEY)
    log.info("T1 residue numbers mapped for %d / %d targets.",
             sum(1 for m in t1_mapping.values()
                 if all(v is not None for v in m.values())),
             len(structures))

    # ── Prepare ligands ────────────────────────────────────────────────────
    lig_pdbqts: dict = {}
    reactive_names: dict = {}
    for lig_name, smiles in LIGANDS.items():
        out = lig_dir / f"{lig_name}.pdbqt"
        if prepare_ligand(lig_name, smiles, out):
            lig_pdbqts[lig_name] = out
            smarts = LIGAND_REACTIVE_SMARTS.get(lig_name)
            if smarts:
                rname = find_reactive_atom_name(out, smarts)
                reactive_names[lig_name] = rname
                log.info("  Reactive atom for %s: %s", lig_name, rname or "NOT FOUND")
        else:
            log.error("Ligand preparation failed for %s – it will be skipped.", lig_name)

    if not lig_pdbqts:
        log.error("No ligands prepared. Aborting.")
        sys.exit(1)

    # ── Positive control (optional) ────────────────────────────────────────
    if args.positive_control and "Guaiacol" in lig_pdbqts:
        if args.dry_run:
            log.info("DRY RUN – would run positive control on %s", ref_pdb.name)
        else:
            run_positive_control(
                ref_pdb, lig_pdbqts["Guaiacol"], prep_dir, res_dir,
                args.box_size, args.exhaustiveness, args.n_poses, args.cpus,
            )

    # ── Main docking loop ──────────────────────────────────────────────────
    summary_records: list = []
    total = len(structures) * len(lig_pdbqts)
    done  = 0

    for uid, pdb_path in structures.items():
        log.info("═" * 60)
        log.info("Protein: %s  (%s)", uid, pdb_path.name)

        # Signal-peptide trim (AlphaFold only)
        is_alphafold = uid not in ITL_IDS
        trimmed_path, n_trimmed = trim_n_terminal(
            pdb_path, ref_pdb, trimmed_dir, is_alphafold=is_alphafold,
        )
        if n_trimmed == 0 and is_alphafold:
            log.info("  No N-terminal residues trimmed for %s.", uid)

        # Prepare receptor PDBQT from (possibly trimmed) PDB
        rec_pdbqt = prep_dir / f"{uid}.pdbqt"
        if not prepare_receptor(trimmed_path, rec_pdbqt):
            log.warning("  Skipping %s – receptor preparation failed.", uid)
            done += len(lig_pdbqts)
            continue

        # Locate T1 site in original structure (MSA mapping → superposition fallback)
        log.info("  Locating T1 site …")
        cu_coords = locate_t1_site(pdb_path, ref_pdb, t1_mapping=t1_mapping.get(uid))
        if cu_coords is None:
            log.warning("  Cannot determine T1 site for %s – skipping.", uid)
            done += len(lig_pdbqts)
            continue

        for lig_name, lig_pdbqt in lig_pdbqts.items():
            done += 1
            tag = f"{uid}_{lig_name}"
            log.info("  [%d/%d] Docking %s …", done, total, tag)

            out_pdbqt = res_dir / f"{tag}_poses.pdbqt"
            out_log   = res_dir / f"{tag}.log"

            if args.dry_run:
                log.info("  DRY RUN – would dock %s into %s at (%.2f, %.2f, %.2f)",
                         lig_name, uid, *cu_coords)
                summary_records.append({
                    "protein_id": uid, "ligand": lig_name,
                    "best_affinity_mean": "DRY_RUN",
                    "best_affinity_sd":   "DRY_RUN",
                    "best_affinity_kcal_mol": "DRY_RUN",
                    "all_affinities_kcal_mol": "DRY_RUN",
                    "t1_cu_x": round(cu_coords[0], 2),
                    "t1_cu_y": round(cu_coords[1], 2),
                    "t1_cu_z": round(cu_coords[2], 2),
                    "reactive_atom_dist_to_T1": "DRY_RUN",
                    "productive_pose": "DRY_RUN",
                    "ranking_score":   "DRY_RUN",
                    "poses_pdbqt": str(out_pdbqt),
                })
                continue

            seed1_results, mean_aff, sd_aff = run_vina_replicates(
                receptor_pdbqt=rec_pdbqt,
                ligand_pdbqt=lig_pdbqt,
                center=cu_coords,
                box_size=args.box_size,
                out_pdbqt=out_pdbqt,
                out_log=out_log,
                n_poses=args.n_poses,
                cpus=args.cpus,
            )

            if seed1_results is None:
                log.error("  Docking failed for %s.", tag)
                summary_records.append({
                    "protein_id": uid, "ligand": lig_name,
                    "best_affinity_mean": "FAILED",
                    "best_affinity_sd":   "FAILED",
                    "best_affinity_kcal_mol": "FAILED",
                    "all_affinities_kcal_mol": "FAILED",
                    "t1_cu_x": round(cu_coords[0], 2),
                    "t1_cu_y": round(cu_coords[1], 2),
                    "t1_cu_z": round(cu_coords[2], 2),
                    "reactive_atom_dist_to_T1": "",
                    "productive_pose": "",
                    "ranking_score": "",
                    "poses_pdbqt": str(out_pdbqt),
                })
                continue

            best_seed1 = min(a for _, a in seed1_results)
            all_aff    = ";".join(f"{a:.3f}" for _, a in sorted(seed1_results))
            log.info("  Seeds mean=%.2f SD=%.2f kcal/mol", mean_aff, sd_aff)

            # T1 geometry filter
            reactive_name = reactive_names.get(lig_name)
            top_dist      = None
            ranking_score = best_seed1
            productive    = False

            if reactive_name and out_pdbqt.exists():
                top_dist = measure_reactive_dist(out_pdbqt, reactive_name, cu_coords, pose_idx=0)

                # Find best affinity among poses with reactive atom ≤ 8 Å
                productive_best: Optional[float] = None
                for pose_idx, (pose_num, aff) in enumerate(sorted(seed1_results, key=lambda x: x[1])):
                    d = measure_reactive_dist(out_pdbqt, reactive_name, cu_coords, pose_idx=pose_idx)
                    if d is not None and d <= 8.0:
                        if productive_best is None or aff < productive_best:
                            productive_best = aff

                if productive_best is not None:
                    ranking_score = productive_best
                    productive    = True
                else:
                    productive = False

                dist_str = f"{top_dist:.2f}" if top_dist is not None else "N/A"
                log.info("  Top pose reactive dist=%.1f Å  productive=%s  ranking_score=%.2f",
                         top_dist if top_dist is not None else float("nan"),
                         productive, ranking_score)

            summary_records.append({
                "protein_id": uid,
                "ligand":     lig_name,
                "best_affinity_mean": round(mean_aff, 3),
                "best_affinity_sd":   round(sd_aff, 3),
                "best_affinity_kcal_mol": round(best_seed1, 3),
                "all_affinities_kcal_mol": all_aff,
                "t1_cu_x": round(cu_coords[0], 2),
                "t1_cu_y": round(cu_coords[1], 2),
                "t1_cu_z": round(cu_coords[2], 2),
                "reactive_atom_dist_to_T1": round(top_dist, 2) if top_dist is not None else "",
                "productive_pose": productive,
                "ranking_score":   round(ranking_score, 3),
                "poses_pdbqt": str(out_pdbqt),
            })

    # ── Write summary ──────────────────────────────────────────────────────
    summary_csv = res_dir / "docking_summary.csv"
    write_summary(summary_records, summary_csv)

    log.info("Done.  %d docking runs recorded.", len(summary_records))
    log.info("Summary CSV: %s", summary_csv)


if __name__ == "__main__":
    main()
