"""
classify_upos.py
================
Step 2: Classify each query as long-UPO or short-UPO using a two-tier approach.

Tier 1 — Sequence-based (fast, no structure needed):
  - Fetch mature protein sequence from UniProt
  - Align to PaDa-I and MroUPO reference sequences
  - Use alignment score + sequence length

Tier 2 — Structural (TM-score, runs after structures are available):
  - Align query structure to both 5OXU and 5FUJ using TMalign
  - Assign based on higher TM-score
  - Flag ambiguous cases (delta TM-score < threshold)

Outputs:
  - results/classification.tsv
  - results/classification_summary.txt
"""

import os
import re
import json
import logging
import subprocess
import tempfile
import requests
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Reference sequences (hardcoded as fallback; fetched from UniProt at runtime)
# ---------------------------------------------------------------------------

REFERENCE_UNIPROT = {
    "long_upo": "B9WCX3",   # AaeUPO (AaeUPO1, basis for PaDa-I)
    "short_upo": "A0A023H437",  # MroUPO (Marasmius rotula) — fragment; fallback below
}

# Full MroUPO sequence from PDB 5FUK (chain A) — used when UniProt fetch returns a fragment
MROUPO_FALLBACK_SEQ = (
    "SAHPWKAPGPNDSRGPCPGLNTLANHGFLPRNGRNISVPMIVKAGFEGYNVQSDILILAGKIGMLTSREA"
    "DTISLEDLKLHGTIEHDASLSREDVAIGDNLHFNEAIFTTLANSNPGADVYNISSAAQVQHDRLADSLAR"
    "NPNVTNTDLTATIRSSESAFFLTVMSAGDPLRGEAPKKFVNVFFREERMPIKEGWKRSTTPITIPLLGPI"
    "IERITELSDWKPTGDNCGAIVLSP"
)

# Approximate mature sequence length ranges (post signal peptide cleavage)
LENGTH_RANGES = {
    "long_upo":  (300, 420),
    "short_upo": (220, 310),
}


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class ClassificationResult:
    accession: str
    seq_length: Optional[int] = None
    length_call: Optional[str] = None      # "long", "short", "ambiguous"
    align_score_long: Optional[float] = None
    align_score_short: Optional[float] = None
    align_identity_long: Optional[float] = None
    align_identity_short: Optional[float] = None
    sequence_call: Optional[str] = None    # "long", "short", "ambiguous"
    tmscore_long: Optional[float] = None
    tmscore_short: Optional[float] = None
    structural_call: Optional[str] = None  # "long", "short", "ambiguous", "no_structure"
    final_call: Optional[str] = None       # consensus
    confidence: str = "low"               # "high", "medium", "low"
    notes: str = ""


# ---------------------------------------------------------------------------
# Sequence fetching
# ---------------------------------------------------------------------------

def fetch_uniprot_sequence(accession: str) -> Optional[str]:
    """Fetch canonical sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        seq = "".join(l for l in lines if not l.startswith(">"))
        return seq
    except Exception as e:
        logger.warning(f"[{accession}] Sequence fetch error: {e}")
        return None


def fetch_reference_sequences(ref_config: dict) -> dict[str, str]:
    """Fetch long and short UPO reference sequences from UniProt."""
    refs = {}
    for key, uniprot_id in REFERENCE_UNIPROT.items():
        seq = fetch_uniprot_sequence(uniprot_id)
        if seq and len(seq) >= 200:
            refs[key] = seq
            logger.info(f"Reference sequence [{key}] {uniprot_id}: {len(seq)} aa")
        else:
            if key == "short_upo":
                refs[key] = MROUPO_FALLBACK_SEQ
                logger.warning(
                    f"UniProt fetch for {key} ({uniprot_id}) returned a fragment or failed — "
                    f"using hardcoded MroUPO fallback ({len(MROUPO_FALLBACK_SEQ)} aa from PDB 5FUK)"
                )
            else:
                logger.error(f"Could not fetch reference sequence for {key} ({uniprot_id})")
    return refs


# ---------------------------------------------------------------------------
# Sequence-level classification
# ---------------------------------------------------------------------------

def classify_by_length(seq_length: int, config: dict) -> str:
    """Quick length-based pre-classification."""
    long_min, long_max = config["classification"]["length_long_min"], \
                         config["classification"]["length_long_max"]
    short_min, short_max = config["classification"]["length_short_min"], \
                           config["classification"]["length_short_max"]

    in_long = long_min <= seq_length <= long_max
    in_short = short_min <= seq_length <= short_max

    if in_long and not in_short:
        return "long"
    elif in_short and not in_long:
        return "short"
    elif in_long and in_short:
        return "ambiguous"  # overlap zone
    else:
        return "ambiguous"  # outside both ranges — unusual, flag


def pairwise_align_score(query_seq: str, ref_seq: str) -> tuple[float, float]:
    """
    Smith-Waterman local alignment.
    Returns (raw_score, percent_identity).
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(query_seq, ref_seq))
    if not alignments:
        return 0.0, 0.0

    best = alignments[0]
    score = best.score

    # Percent identity: matches / aligned positions
    aln_str = str(best).splitlines()  # query / match / target lines
    if len(aln_str) >= 3:
        q_line = aln_str[0]
        t_line = aln_str[2]
        matches = sum(a == b for a, b in zip(q_line, t_line)
                      if a not in ("-", " ") and b not in ("-", " "))
        aligned = sum(1 for a, b in zip(q_line, t_line)
                      if a not in ("-", " ") and b not in ("-", " "))
        identity = (matches / aligned * 100) if aligned > 0 else 0.0
    else:
        identity = 0.0

    return score, identity


def classify_by_sequence(accession: str, query_seq: str,
                         ref_seqs: dict, config: dict) -> ClassificationResult:
    """Full sequence-level classification."""
    result = ClassificationResult(accession=accession)
    result.seq_length = len(query_seq)

    # Length call
    result.length_call = classify_by_length(result.seq_length, config)
    logger.info(f"[{accession}] Length={result.seq_length} aa → {result.length_call}")

    # Alignment scores vs both references
    if "long_upo" in ref_seqs:
        score_l, id_l = pairwise_align_score(query_seq, ref_seqs["long_upo"])
        result.align_score_long = score_l
        result.align_identity_long = id_l

    if "short_upo" in ref_seqs:
        score_s, id_s = pairwise_align_score(query_seq, ref_seqs["short_upo"])
        result.align_score_short = score_s
        result.align_identity_short = id_s

    # Sequence call based on alignment
    if result.align_score_long is not None and result.align_score_short is not None:
        diff = result.align_score_long - result.align_score_short
        if diff > 50:
            result.sequence_call = "long"
        elif diff < -50:
            result.sequence_call = "short"
        else:
            result.sequence_call = "ambiguous"
    elif result.align_score_long is not None:
        result.sequence_call = "long"
    elif result.align_score_short is not None:
        result.sequence_call = "short"
    else:
        result.sequence_call = "ambiguous"

    score_l_str = f"{result.align_score_long:.0f}" if result.align_score_long is not None else "N/A"
    score_s_str = f"{result.align_score_short:.0f}" if result.align_score_short is not None else "N/A"
    logger.info(
        f"[{accession}] Align scores — long:{score_l_str} "
        f"short:{score_s_str} → seq_call:{result.sequence_call}"
    )

    # Partial consensus (sequence tier only)
    if result.length_call == result.sequence_call and result.length_call != "ambiguous":
        result.final_call = result.length_call
        result.confidence = "medium"
    elif result.sequence_call != "ambiguous":
        result.final_call = result.sequence_call
        result.confidence = "medium"
    else:
        result.final_call = result.length_call if result.length_call != "ambiguous" else "ambiguous"
        result.confidence = "low"

    return result


# ---------------------------------------------------------------------------
# Structural classification via TM-align
# ---------------------------------------------------------------------------

def run_tmalign(query_pdb: Path, ref_pdb: Path) -> Optional[float]:
    """
    Run TMalign and extract TM-score (normalised to query length).
    Returns None if TMalign is not available or fails.
    """
    try:
        result = subprocess.run(
            ["TMalign", str(query_pdb), str(ref_pdb)],
            capture_output=True, text=True, timeout=60
        )
        # Parse TM-score (line: "TM-score= X (if normalized by length of Chain_1)")
        for line in result.stdout.splitlines():
            if "TM-score=" in line and "Chain_1" in line:
                match = re.search(r"TM-score=\s*([\d.]+)", line)
                if match:
                    return float(match.group(1))
    except FileNotFoundError:
        logger.warning("TMalign not found in PATH. Skipping structural classification.")
    except subprocess.TimeoutExpired:
        logger.warning(f"TMalign timed out for {query_pdb.name}")
    return None


def classify_by_structure(result: ClassificationResult,
                          structure_path: Path,
                          ref_structures: dict[str, Path],
                          config: dict) -> ClassificationResult:
    """
    Augment an existing ClassificationResult with TM-score based structural call.
    Requires TMalign in PATH.
    """
    if structure_path is None or not structure_path.exists():
        result.structural_call = "no_structure"
        return result

    tmscore_long, tmscore_short = None, None

    if "long_upo" in ref_structures:
        tmscore_long = run_tmalign(structure_path, ref_structures["long_upo"])
        result.tmscore_long = tmscore_long

    if "short_upo" in ref_structures:
        tmscore_short = run_tmalign(structure_path, ref_structures["short_upo"])
        result.tmscore_short = tmscore_short

    # Structural call
    min_diff = config["classification"]["tmscore_min_diff"]
    if tmscore_long is not None and tmscore_short is not None:
        diff = tmscore_long - tmscore_short
        if diff > min_diff:
            result.structural_call = "long"
        elif diff < -min_diff:
            result.structural_call = "short"
        else:
            result.structural_call = "ambiguous"
            result.notes += f" TM-score diff={diff:.3f} below threshold."
    elif tmscore_long is not None:
        result.structural_call = "long" if tmscore_long > 0.5 else "ambiguous"
    elif tmscore_short is not None:
        result.structural_call = "short" if tmscore_short > 0.5 else "ambiguous"
    else:
        result.structural_call = "no_tmalign"

    logger.info(
        f"[{result.accession}] TM-scores — long:{tmscore_long} "
        f"short:{tmscore_short} → struct_call:{result.structural_call}"
    )

    # Update final call with structural evidence
    _update_final_call(result, config)
    return result


def _update_final_call(result: ClassificationResult, config: dict):
    """Consensus logic across all tiers."""
    calls = [c for c in [result.length_call, result.sequence_call,
                          result.structural_call]
             if c and c not in ("ambiguous", "no_structure", "no_tmalign")]

    if not calls:
        result.final_call = "ambiguous"
        result.confidence = "low"
        return

    long_votes = calls.count("long")
    short_votes = calls.count("short")
    total = len(calls)

    if long_votes == total:
        result.final_call = "long"
        result.confidence = "high" if total >= 2 else "medium"
    elif short_votes == total:
        result.final_call = "short"
        result.confidence = "high" if total >= 2 else "medium"
    elif long_votes > short_votes:
        result.final_call = "long"
        result.confidence = "medium"
        result.notes += " Majority vote (long)."
    elif short_votes > long_votes:
        result.final_call = "short"
        result.confidence = "medium"
        result.notes += " Majority vote (short)."
    else:
        result.final_call = "ambiguous"
        result.confidence = "low"
        result.notes += " Tied votes — manual inspection needed."


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_classification(results: list[ClassificationResult], output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / "classification.tsv"

    with open(tsv_path, "w") as f:
        f.write(
            "accession\tseq_length\tlength_call\talign_score_long\talign_score_short\t"
            "align_id_long\talign_id_short\tsequence_call\t"
            "tmscore_long\ttmscore_short\tstructural_call\t"
            "final_call\tconfidence\tnotes\n"
        )
        for r in results:
            f.write(
                f"{r.accession}\t{r.seq_length or ''}\t{r.length_call or ''}\t"
                f"{r.align_score_long or ''}\t{r.align_score_short or ''}\t"
                f"{r.align_identity_long or ''}\t{r.align_identity_short or ''}\t"
                f"{r.sequence_call or ''}\t"
                f"{r.tmscore_long or ''}\t{r.tmscore_short or ''}\t"
                f"{r.structural_call or ''}\t"
                f"{r.final_call or ''}\t{r.confidence}\t{r.notes}\n"
            )

    # Summary
    summary_path = output_dir / "classification_summary.txt"
    long_hi  = [r for r in results if r.final_call == "long"  and r.confidence == "high"]
    long_med = [r for r in results if r.final_call == "long"  and r.confidence == "medium"]
    short_hi = [r for r in results if r.final_call == "short" and r.confidence == "high"]
    short_med= [r for r in results if r.final_call == "short" and r.confidence == "medium"]
    ambig    = [r for r in results if r.final_call == "ambiguous"]

    lines = [
        "UPO CLASSIFICATION SUMMARY",
        "=" * 40,
        f"Total queries: {len(results)}",
        f"Long-UPO  (high confidence): {len(long_hi)}",
        f"Long-UPO  (medium confidence): {len(long_med)}",
        f"Short-UPO (high confidence): {len(short_hi)}",
        f"Short-UPO (medium confidence): {len(short_med)}",
        f"Ambiguous: {len(ambig)}",
        "",
        "LONG-UPO:",
        *[f"  {r.accession} ({r.confidence})" for r in long_hi + long_med],
        "",
        "SHORT-UPO:",
        *[f"  {r.accession} ({r.confidence})" for r in short_hi + short_med],
        "",
        "AMBIGUOUS (manual review):",
        *[f"  {r.accession}: {r.notes}" for r in ambig],
    ]
    summary_path.write_text("\n".join(lines))
    logger.info(f"Classification summary → {summary_path}")

    return tsv_path, summary_path
