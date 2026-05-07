"""
retrieve_structures.py
======================
Step 1: For each UniProt accession, attempt structure retrieval in order:
  1. Experimental structure via PDB (UniProt → PDB cross-references)
  2. AlphaFold DB predicted structure
  3. Flag for Boltz-2 local prediction if both fail

Outputs:
  - Downloaded .pdb / .cif files in structure_dir/
  - results/structure_sources.tsv  (accession, source, pdb_id, resolution, path)
  - results/boltz_todo.txt          (accessions needing Boltz-2 prediction)
"""

import os
import json
import time
import logging
import requests
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class StructureResult:
    accession: str
    source: str                     # "PDB", "AlphaFold", "BOLTZ_NEEDED"
    pdb_id: Optional[str] = None
    chain: Optional[str] = None
    resolution: Optional[float] = None
    method: Optional[str] = None    # X-RAY, EM, NMR, PREDICTED
    file_path: Optional[Path] = None
    mean_plddt: Optional[float] = None  # AlphaFold only
    notes: str = ""


# ---------------------------------------------------------------------------
# PDB retrieval
# ---------------------------------------------------------------------------

def get_pdb_ids_for_uniprot(accession: str) -> list[dict]:
    """
    Query UniProt → PDB cross-references via the UniProt REST API.
    Returns list of dicts with keys: pdb_id, chain, resolution, method.
    Sorted by resolution ascending (best first), NMR/EM deprioritised.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"[{accession}] UniProt API error: {e}")
        return []

    data = resp.json()
    pdb_refs = []

    for db_ref in data.get("uniProtKBCrossReferences", []):
        if db_ref.get("database") != "PDB":
            continue
        pdb_id = db_ref.get("id")
        method, resolution, chains = None, None, None
        for prop in db_ref.get("properties", []):
            if prop["key"] == "Method":
                method = prop["value"]
            elif prop["key"] == "Resolution":
                try:
                    resolution = float(prop["value"].replace(" A", "").strip())
                except ValueError:
                    pass
            elif prop["key"] == "Chains":
                chains = prop["value"]

        pdb_refs.append({
            "pdb_id": pdb_id,
            "chain": chains,
            "resolution": resolution,
            "method": method,
        })

    # Sort: prefer X-RAY with best resolution
    def sort_key(r):
        method_rank = 0 if r["method"] == "X-ray" else (1 if r["method"] == "EM" else 2)
        res = r["resolution"] if r["resolution"] is not None else 99.0
        return (method_rank, res)

    pdb_refs.sort(key=sort_key)
    return pdb_refs


def download_pdb(pdb_id: str, out_dir: Path, prefer_cif: bool = False) -> Optional[Path]:
    """Download structure from RCSB. Tries mmCIF then PDB format."""
    pdb_id_lower = pdb_id.lower()
    out_dir.mkdir(parents=True, exist_ok=True)

    formats = [("cif", "mmCIF")] if prefer_cif else [("pdb", "PDB"), ("cif", "mmCIF")]
    for fmt, label in formats:
        out_path = out_dir / f"{pdb_id_lower}.{fmt}"
        if out_path.exists():
            logger.info(f"  [{pdb_id}] Already downloaded ({fmt})")
            return out_path
        url = f"https://files.rcsb.org/download/{pdb_id_lower}.{fmt}"
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                out_path.write_bytes(resp.content)
                logger.info(f"  [{pdb_id}] Downloaded {label} → {out_path.name}")
                return out_path
        except requests.RequestException as e:
            logger.warning(f"  [{pdb_id}] Download error ({label}): {e}")

    return None


def retrieve_from_pdb(accession: str, out_dir: Path,
                      max_resolution: float = 3.0) -> Optional[StructureResult]:
    """
    Full PDB retrieval flow for one accession.
    Skips structures worse than max_resolution Å (still reports NMR/EM).
    """
    pdb_refs = get_pdb_ids_for_uniprot(accession)
    if not pdb_refs:
        logger.info(f"[{accession}] No PDB cross-references found.")
        return None

    for ref in pdb_refs:
        pdb_id = ref["pdb_id"]
        res = ref["resolution"]
        method = ref["method"]

        # Skip poor-resolution X-ray structures but keep NMR/EM
        if method == "X-ray" and res is not None and res > max_resolution:
            logger.info(f"  [{accession}] Skipping {pdb_id} (res={res} Å > {max_resolution} Å)")
            continue

        logger.info(f"  [{accession}] Trying PDB:{pdb_id} method={method} res={res}")
        file_path = download_pdb(pdb_id, out_dir)
        if file_path:
            return StructureResult(
                accession=accession,
                source="PDB",
                pdb_id=pdb_id,
                chain=ref["chain"],
                resolution=res,
                method=method,
                file_path=file_path,
                notes=f"Best PDB entry ({len(pdb_refs)} total cross-refs)"
            )

    logger.info(f"[{accession}] All PDB entries failed or exceeded resolution cutoff.")
    return None


# ---------------------------------------------------------------------------
# AlphaFold DB retrieval
# ---------------------------------------------------------------------------

def retrieve_from_alphafold(accession: str, out_dir: Path,
                            version: int = 4) -> Optional[StructureResult]:
    """
    Download AlphaFold DB structure for a UniProt accession.
    Returns StructureResult or None if not found.
    Also extracts mean pLDDT from the B-factor column as a QC metric.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"AF-{accession}-F1-model_v{version}.pdb"

    if not out_path.exists():
        url = (f"https://alphafold.ebi.ac.uk/files/"
               f"AF-{accession}-F1-model_v{version}.pdb")
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 404:
                logger.info(f"[{accession}] Not in AlphaFold DB (v{version}).")
                return None
            resp.raise_for_status()
            out_path.write_bytes(resp.content)
            logger.info(f"[{accession}] AlphaFold DB → {out_path.name}")
        except requests.RequestException as e:
            logger.warning(f"[{accession}] AlphaFold download error: {e}")
            return None
    else:
        logger.info(f"[{accession}] AlphaFold already downloaded.")

    # Extract mean pLDDT from B-factor column
    mean_plddt = _extract_mean_plddt(out_path)

    return StructureResult(
        accession=accession,
        source="AlphaFold",
        pdb_id=f"AF-{accession}-F1-v{version}",
        chain="A",
        resolution=None,
        method="PREDICTED",
        file_path=out_path,
        mean_plddt=mean_plddt,
        notes=f"AlphaFold DB v{version}, mean pLDDT={mean_plddt:.1f}" if mean_plddt else ""
    )


def _extract_mean_plddt(pdb_path: Path) -> Optional[float]:
    """Parse B-factor column of CA atoms as pLDDT proxy."""
    try:
        bfactors = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        bfactors.append(float(line[60:66].strip()))
                    except ValueError:
                        pass
        if bfactors:
            return sum(bfactors) / len(bfactors)
    except Exception as e:
        logger.warning(f"pLDDT extraction failed for {pdb_path.name}: {e}")
    return None


# ---------------------------------------------------------------------------
# Reference structure download
# ---------------------------------------------------------------------------

def download_references(config: dict, out_dir: Path) -> dict[str, Path]:
    """Download PaDa-I (5OXU) and MroUPO (5FUJ) reference structures."""
    refs = {}
    for key, ref_cfg in config["references"].items():
        pdb_id = ref_cfg["pdb_id"]
        path = download_pdb(pdb_id, out_dir / "references")
        if path:
            refs[key] = path
            logger.info(f"Reference [{key}] {pdb_id} → {path.name}")
        else:
            logger.error(f"FAILED to download reference structure {pdb_id}")
    return refs


# ---------------------------------------------------------------------------
# Main retrieval loop
# ---------------------------------------------------------------------------

def run_retrieval(config: dict) -> list[StructureResult]:
    """
    Run full retrieval pipeline for all accessions in config.
    Returns list of StructureResult (one per accession).
    """
    structure_dir = Path(config["structure_dir"])
    structure_dir.mkdir(parents=True, exist_ok=True)

    results = []
    boltz_todo = []

    for acc in config["queries"]:
        logger.info(f"\n{'='*50}")
        logger.info(f"Processing: {acc}")

        # 1. Try PDB
        result = retrieve_from_pdb(acc, structure_dir / "pdb")
        if result:
            results.append(result)
            time.sleep(0.5)  # be polite to RCSB
            continue

        # 2. Try AlphaFold DB
        time.sleep(0.3)
        result = retrieve_from_alphafold(acc, structure_dir / "alphafold")
        if result:
            results.append(result)
            # Warn if predicted structure has low pLDDT
            if result.mean_plddt and result.mean_plddt < 70:
                logger.warning(
                    f"[{acc}] Low mean pLDDT={result.mean_plddt:.1f} — "
                    "consider Boltz-2 re-prediction or manual inspection"
                )
            continue

        # 3. Flag for Boltz-2
        logger.warning(f"[{acc}] No structure found — flagged for Boltz-2")
        boltz_todo.append(acc)
        results.append(StructureResult(
            accession=acc,
            source="BOLTZ_NEEDED",
            notes="Not found in PDB or AlphaFold DB"
        ))

    return results, boltz_todo


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_results(results: list[StructureResult], boltz_todo: list[str],
                  output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    # TSV summary
    tsv_path = output_dir / "structure_sources.tsv"
    with open(tsv_path, "w") as f:
        f.write("accession\tsource\tpdb_id\tchain\tresolution\tmethod\t"
                "mean_plddt\tfile_path\tnotes\n")
        for r in results:
            f.write(
                f"{r.accession}\t{r.source}\t{r.pdb_id or ''}\t"
                f"{r.chain or ''}\t{r.resolution or ''}\t{r.method or ''}\t"
                f"{r.mean_plddt or ''}\t{r.file_path or ''}\t{r.notes}\n"
            )
    logger.info(f"Structure sources → {tsv_path}")

    # Boltz todo list
    if boltz_todo:
        boltz_path = output_dir / "boltz_todo.txt"
        boltz_path.write_text("\n".join(boltz_todo) + "\n")
        logger.info(f"Boltz-2 queue ({len(boltz_todo)} entries) → {boltz_path}")
        print(f"\n⚠  {len(boltz_todo)} accession(s) need Boltz-2 prediction:")
        for acc in boltz_todo:
            print(f"   {acc}")
        print(f"   Run: boltz predict <fasta> --out_dir results/structures/boltz/")

    return tsv_path
