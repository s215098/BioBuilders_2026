#!/usr/bin/env python3
"""
graft_heme.py
=============
Graft heme from a reference UPO crystal structure onto predicted (Boltz-2 /
AlphaFold) structures by structural superposition.

Strategy
--------
1. Superimpose the predicted structure (CA atoms) onto the appropriate reference
   (5OXU for long-UPO, 5FUJ chain A for short-UPO).
2. Apply the same rotation/translation to the reference heme atoms.
3. Insert the transformed heme as a new residue in the predicted structure.
4. Verify that the Cys axial ligand S–Fe distance is ≤ 3.5 Å.
5. Write the grafted structure to results/structures/grafted/.

Usage
-----
    python graft_heme.py --config config.yaml
    python graft_heme.py --config config.yaml --accession A0A8T9CJ62
"""

import argparse
import csv
import logging
import sys
import numpy as np
from copy import deepcopy
from pathlib import Path

import yaml
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
from Bio.PDB import Selection, PPBuilder
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)
PARSER = PDBParser(QUIET=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_model(path: Path, sid: str):
    struct = PARSER.get_structure(sid, str(path))
    return struct[0]


AA3TO1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "MSE":"M","SEP":"S","TPO":"T","CSO":"C",
}


def get_ca_atoms(model, chain_id: str = None) -> list:
    """Return all CA atoms (with residue), optionally from a single chain."""
    atoms = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            if "CA" in res:
                atoms.append(res["CA"])
    return atoms


def get_residues_with_ca(model, chain_id: str = None) -> list:
    """Return list of (residue, CA_atom) for standard residues."""
    pairs = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            if "CA" in res:
                pairs.append((res, res["CA"]))
    return pairs


def sequence_from_residues(res_ca_pairs: list) -> str:
    return "".join(AA3TO1.get(r.resname, "X") for r, _ in res_ca_pairs)


def matched_ca_pairs(query_pairs: list, ref_pairs: list) -> tuple[list, list]:
    """
    Align query and ref sequences; return matched (query_CA, ref_CA) lists.
    Only aligned (non-gap) positions are included.
    """
    q_seq = sequence_from_residues(query_pairs)
    r_seq = sequence_from_residues(ref_pairs)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score  = -10
    aligner.extend_gap_score = -0.5

    alignment = next(iter(aligner.align(q_seq, r_seq)))
    q_indices, r_indices = [], []
    qi, ri = 0, 0
    for q_char, r_char in zip(*alignment):
        if q_char != "-" and r_char != "-":
            q_indices.append(qi)
            r_indices.append(ri)
        if q_char != "-":
            qi += 1
        if r_char != "-":
            ri += 1

    q_cas = [query_pairs[i][1] for i in q_indices]
    r_cas = [ref_pairs[i][1]   for i in r_indices]
    return q_cas, r_cas


def get_axial_cys(model, chain_id: str = None, fe_coord: np.ndarray = None,
                  max_dist: float = 5.0):
    """
    Return the Cys residue closest to fe_coord (axial ligand).
    If fe_coord is None, return the first Cys in the structure.
    """
    best_res, best_dist = None, max_dist
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for res in chain:
            if res.resname != "CYS":
                continue
            if fe_coord is not None and "SG" in res:
                d = np.linalg.norm(res["SG"].coord - fe_coord)
                if d < best_dist:
                    best_res, best_dist = res, d
            elif fe_coord is None and best_res is None:
                best_res = res
    return best_res


def cys_anchor_atoms(cys_res) -> list:
    """
    Return [N, CA, C, CB, SG] atoms of a Cys residue (skip missing atoms).
    These 5 atoms define a rigid local frame for heme placement.
    """
    names = ["N", "CA", "C", "CB", "SG"]
    return [cys_res[n] for n in names if n in cys_res]


def get_heme_residue(model, chain_id: str = "A"):
    """Return the first HEM residue in the given chain."""
    for chain in model:
        if chain.id != chain_id:
            continue
        for res in chain:
            if res.resname in ("HEM", "HEC", "HEO"):
                return res
    return None


def find_cys_axial(model, fe_coord: np.ndarray, max_dist: float = 3.5):
    """Return (residue, S-atom) of the Cys whose SG is closest to Fe."""
    best, best_dist = None, max_dist
    for chain in model:
        for res in chain:
            if res.resname != "CYS":
                continue
            for atom in res:
                if atom.name in ("SG", "S"):
                    d = np.linalg.norm(atom.coord - fe_coord)
                    if d < best_dist:
                        best = (res, atom)
                        best_dist = d
    return best, best_dist


def transform_heme(heme_res, rot: np.ndarray, tran: np.ndarray):
    """Return a deep copy of heme_res with all atom coords transformed."""
    heme_copy = deepcopy(heme_res)
    for atom in heme_copy:
        atom.coord = rot @ atom.coord + tran
    return heme_copy


class HemeSelect(Select):
    """PDBIO selector that writes all ATOM records plus a single heme residue."""
    def __init__(self, heme_res):
        self.heme_id = (heme_res.get_parent().id, heme_res.id)

    def accept_residue(self, res):
        chain_id = res.get_parent().id
        if res.id[0] == " ":          # standard amino acid
            return True
        if (chain_id, res.id) == self.heme_id:
            return True
        return False


# ---------------------------------------------------------------------------
# Core grafting
# ---------------------------------------------------------------------------

def graft_heme(
    query_path: Path,
    ref_path: Path,
    ref_chain: str,
    out_path: Path,
    min_ca_overlap: int = 3,
) -> dict:
    """
    Superimpose query onto ref, transfer heme, write grafted PDB.

    Returns a result dict with keys:
        success, ca_rmsd, fe_coord, cys_res, cys_dist, flags
    """
    result = {"success": False, "flags": []}

    query_model = load_model(query_path, query_path.stem)
    ref_model   = load_model(ref_path,   ref_path.stem)

    # --- locate heme and axial Cys in reference ----------------------------
    ref_heme_pre = get_heme_residue(ref_model, chain_id=ref_chain)
    if ref_heme_pre is None:
        result["flags"].append("NO_REF_HEME")
        logger.error("No HEM residue found in reference structure")
        return result

    ref_fe_coord = None
    for atom in ref_heme_pre:
        if atom.element == "FE" or atom.name in ("FE", "FE2"):
            ref_fe_coord = np.array(atom.coord)
            break

    ref_cys = get_axial_cys(ref_model, chain_id=ref_chain, fe_coord=ref_fe_coord)
    if ref_cys is None:
        result["flags"].append("NO_REF_CYS_AXIAL")
        logger.error("No axial Cys found in reference near Fe")
        return result

    query_cys = get_axial_cys(query_model)
    if query_cys is None:
        result["flags"].append("NO_QUERY_CYS")
        logger.error("No Cys found in query structure")
        return result

    logger.info(f"Ref axial Cys: {ref_cys.id[1]} chain {ref_cys.get_parent().id}")
    logger.info(f"Query axial Cys: {query_cys.id[1]} chain {query_cys.get_parent().id}")

    # --- Cys-anchored superposition ----------------------------------------
    # Superimpose backbone+SG of ref axial Cys onto query axial Cys.
    # query = fixed, ref = moving → rotran brings ref heme into query space.
    ref_anchor   = cys_anchor_atoms(ref_cys)
    query_anchor = cys_anchor_atoms(query_cys)
    n = min(len(ref_anchor), len(query_anchor))

    if n < min_ca_overlap:
        result["flags"].append(f"TOO_FEW_ANCHOR({n})")
        logger.warning(f"Only {n} anchor atoms — skipping")
        return result

    sup = Superimposer()
    sup.set_atoms(query_anchor[:n], ref_anchor[:n])  # query fixed, ref moving
    ca_rmsd = round(sup.rms, 3)
    result["ca_rmsd"] = ca_rmsd
    logger.info(f"Cys-anchor RMSD: {ca_rmsd} Å  (n={n} atoms)")

    # --- transform heme into query space via sup.apply() -------------------
    grafted_heme = deepcopy(ref_heme_pre)
    sup.apply(list(grafted_heme.get_atoms()))

    # Insert into query model's first chain (predicted structures are single-chain)
    first_chain = next(iter(query_model))
    # Avoid residue ID collision
    existing_ids = {r.id for r in first_chain}
    heme_id = ("H_HEM", 999, " ")
    while heme_id in existing_ids:
        heme_id = ("H_HEM", heme_id[1] + 1, " ")
    grafted_heme.id = heme_id
    grafted_heme.detach_parent()
    first_chain.add(grafted_heme)

    # --- verify Cys axial ligand -------------------------------------------
    fe_coord = None
    for atom in grafted_heme:
        if atom.element == "FE" or atom.name in ("FE", "FE2"):
            fe_coord = np.array(atom.coord)
            break

    if fe_coord is None:
        result["flags"].append("NO_FE_IN_HEME")
        logger.warning("No FE atom found in grafted heme")
    else:
        result["fe_coord"] = tuple(float(x) for x in fe_coord)
        cys_hit, cys_dist = find_cys_axial(query_model, fe_coord)
        result["cys_dist"] = round(cys_dist, 3)
        if cys_hit:
            res, atom = cys_hit
            result["cys_res"] = f"CYS_{res.id[1]}_{res.get_parent().id}"
            if cys_dist <= 3.5:
                logger.info(f"Cys axial ligand: {result['cys_res']}  S–Fe = {cys_dist:.2f} Å ✓")
            else:
                result["flags"].append(f"CYS_FE_DIST_HIGH({cys_dist:.2f})")
                logger.warning(f"S–Fe distance {cys_dist:.2f} Å > 3.5 Å — check alignment")
        else:
            result["flags"].append("NO_CYS_AXIAL")
            logger.warning("No Cys found near grafted Fe")

    # --- write output ------------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(query_model)
    io.save(str(out_path), select=HemeSelect(grafted_heme))
    logger.info(f"Grafted structure → {out_path}")

    result["success"] = True
    return result


# ---------------------------------------------------------------------------
# Pipeline driver
# ---------------------------------------------------------------------------

def run_grafting(config: dict):
    output_dir  = Path(config["output_dir"])
    struct_dir  = Path(config["structure_dir"])
    grafted_dir = struct_dir / "grafted"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    ref_paths = {
        "long_upo":  struct_dir / "references" / "5oxu.pdb",
        "short_upo": struct_dir / "references" / "5fuj.pdb",
    }
    ref_chains = {"long_upo": "A", "short_upo": "A"}

    # Load classification
    clf_map = {}
    clf_tsv = output_dir / "classification.tsv"
    if clf_tsv.exists():
        with open(clf_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                clf_map[row["accession"]] = row["final_call"]

    # Load structure sources — only predicted structures need grafting
    sources_tsv = output_dir / "structure_sources.tsv"
    if not sources_tsv.exists():
        print("ERROR: Run step 1 first (structure_sources.tsv not found)")
        sys.exit(1)

    results = []
    with open(sources_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc    = row["accession"]
            source = row["source"]
            fpath  = row.get("file_path", "").strip()

            if source in ("PDB",):
                print(f"  [{acc}] PDB structure — no heme grafting needed")
                continue
            if not fpath or not Path(fpath).exists():
                print(f"  [{acc}] No structure available — skipping")
                continue

            call = clf_map.get(acc, "long")
            ref_key = "short_upo" if call == "short" else "long_upo"
            ref_path = ref_paths[ref_key]

            if not ref_path.exists():
                print(f"  [{acc}] Reference {ref_path.name} not found — skipping")
                continue

            print(f"\n  [{acc}] source={source} call={call} ref={ref_path.name}")
            out_path = grafted_dir / f"{acc}_grafted.pdb"

            res = graft_heme(
                query_path=Path(fpath),
                ref_path=ref_path,
                ref_chain=ref_chains[ref_key],
                out_path=out_path,
            )
            res["accession"] = acc
            res["source"]    = source
            res["ref"]       = ref_path.name
            results.append(res)

    # Write summary TSV
    tsv_path = output_dir / "heme_grafting.tsv"
    with open(tsv_path, "w") as f:
        f.write("accession\tsource\tref\tsuccess\tca_rmsd\tfe_coord\tcys_res\tcys_dist\tflags\n")
        for r in results:
            f.write(
                f"{r['accession']}\t{r['source']}\t{r.get('ref','')}\t{r['success']}\t"
                f"{r.get('ca_rmsd','')}\t{r.get('fe_coord','')}\t"
                f"{r.get('cys_res','')}\t{r.get('cys_dist','')}\t"
                f"{';'.join(r.get('flags',[]))}\n"
            )
    print(f"\n  Summary → {tsv_path}")

    passed = [r for r in results if r["success"] and not r.get("flags")]
    flagged = [r for r in results if r["success"] and r.get("flags")]
    failed  = [r for r in results if not r["success"]]
    print(f"  Clean: {len(passed)} | Flagged: {len(flagged)} | Failed: {len(failed)}")
    for r in flagged + failed:
        print(f"  {'✗' if not r['success'] else '⚠'} {r['accession']}: {r.get('flags')}")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Graft heme onto predicted UPO structures")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--accession", help="Run for a single accession only")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    with open(args.config) as f:
        config = yaml.safe_load(f)

    print("\n" + "="*60)
    print("HEME GRAFTING")
    print("="*60)

    if args.accession:
        # single-accession mode
        struct_dir = Path(config["structure_dir"])
        output_dir = Path(config["output_dir"])
        clf_map = {}
        clf_tsv = output_dir / "classification.tsv"
        if clf_tsv.exists():
            with open(clf_tsv) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    clf_map[row["accession"]] = row["final_call"]

        sources_tsv = output_dir / "structure_sources.tsv"
        fpath = None
        source = None
        with open(sources_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row["accession"] == args.accession:
                    fpath  = row.get("file_path", "").strip()
                    source = row["source"]
                    break

        if not fpath or not Path(fpath).exists():
            print(f"No structure found for {args.accession}")
            sys.exit(1)

        call    = clf_map.get(args.accession, "long")
        ref_key = "short_upo" if call == "short" else "long_upo"
        ref_path = struct_dir / "references" / ("5fuj.pdb" if ref_key == "short_upo" else "5oxu.pdb")
        out_path = struct_dir / "grafted" / f"{args.accession}_grafted.pdb"

        res = graft_heme(Path(fpath), ref_path, "A", out_path)
        print(res)
    else:
        run_grafting(config)

    print("\n✓ Done.")


if __name__ == "__main__":
    main()
