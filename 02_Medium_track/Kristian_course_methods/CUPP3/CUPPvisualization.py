#!/usr/bin/env python3
"""
CUPPvisualization.py

Usage example:
  python3 CUPPvisualization.py \
      -pool CUPP/CUPPpools/GH29_CUPPpool.json \
      -tree CUPP/itol/GH29_fa8x2_90.tree \
      -name GH29 \
      -outdir CUPP/itol

Requirements:
  pip install biopython
"""

import os
import json
import argparse
import random
from typing import Union, List, Dict, Tuple, Set

# Use BioPython for safe Newick parsing and pruning
try:
    from Bio import Phylo
except ImportError as e:
    raise SystemExit(
        "BioPython is required. Install with: pip install biopython\n"
        f"Original import error: {e}"
    )


# ------------------------------
# Helpers
# ------------------------------

import requests, json, os

def load_accessions_from_fasta(path):
    """
    Extract accessions from FASTA headers supporting multiple formats:
    1. UniProt format: >tr|A0ABX5ZX28|A0ABX5ZX28_STRTE Hyaluronidase...
    2. Simple format: >D3GBU0 Putative laccase 5 OS=Agaricus...
    3. Other formats with pipes: >prefix|ACCESSION|suffix
    """
    accessions = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = None
                
                # Check if it's a pipe-delimited format
                if "|" in header:
                    parts = header.split("|")
                    if len(parts) >= 2:
                        # For UniProt format like tr|A0ABX5ZX28|...
                        # Take the second part as accession
                        acc = parts[1].split()[0]
                    else:
                        # Single pipe case, take first part
                        acc = parts[0].split()[0]
                else:
                    # Simple format: take first word as accession
                    acc = header.split()[0]
                
                if acc:
                    accessions.append(acc)
    return accessions


def fetch_uniprot_metadata(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None
    return r.json()


def extract_pfam(metadata):
    pfam_entries = []
    if metadata is None:
        return pfam_entries
    xrefs = metadata.get("uniProtKBCrossReferences", [])
    for x in xrefs:
        if x.get("database") == "Pfam":
            pfam_id = x.get("id")
            pfam_entries.append(pfam_id)
    return pfam_entries


# -----------------------------
# NEW FUNCTION: Only download if missing
# -----------------------------
def load_or_fetch_metadata(accessions, out_json):
    # Load existing metadata if available
    existing_metadata = {}
    if os.path.exists(out_json) and os.path.getsize(out_json) > 0:
        try:
            with open(out_json, "r") as f:
                print(f"Loading existing metadata from {out_json}")
                existing_metadata = json.load(f)
            print(f"✓ Loaded {len(existing_metadata)} existing metadata entries")
        except json.JSONDecodeError:
            print(f"WARNING: Corrupted metadata file {out_json}, starting fresh")
            existing_metadata = {}
        except Exception as e:
            print(f"ERROR loading metadata: {e}")
            existing_metadata = {}
    else:
        print(f"No existing metadata found, will create {out_json}")
    
    # Check which accessions are missing from existing metadata
    missing_accessions = []
    for acc in accessions:
        # Check both acc and acc_SPECIES formats (e.g., O60502 and O60502_HUMAN)
        found = False
        for key in existing_metadata.keys():
            if key.startswith(acc) or key.split("_")[0] == acc:
                found = True
                break
        if not found:
            missing_accessions.append(acc)
    
    print(f"Accession analysis:")
    print(f"  Total accessions from FASTA: {len(accessions)}")
    print(f"  Already in metadata: {len(accessions) - len(missing_accessions)}")
    print(f"  Missing from metadata: {len(missing_accessions)}")
    
    # If all accessions are covered, return existing metadata
    if not missing_accessions:
        print(f"✓ All {len(accessions)} proteins found in existing metadata")
        return existing_metadata
    
    # Download metadata for missing accessions
    print(f"Fetching {len(missing_accessions)} missing proteins from UniProt...")
    if len(missing_accessions) > 20:
        print(f"  This may take a while for {len(missing_accessions)} proteins...")
    
    results = existing_metadata.copy()  # Start with existing data
    successful_fetches = 0
    failed_fetches = 0
    
    for i, acc in enumerate(missing_accessions, 1):
        print(f"Fetching {acc}... ({i}/{len(missing_accessions)})")
        try:
            meta = fetch_uniprot_metadata(acc)
            if meta:
                pfam = extract_pfam(meta)
                # Use UniProt ID as key (e.g., O60502_HUMAN)
                uniprot_id = meta.get("uniProtkbId", acc)
                results[uniprot_id] = {
                    "pfam": pfam,
                    "raw": meta
                }
                print(f"✓ {acc} -> {uniprot_id}")
                successful_fetches += 1
            else:
                print(f"✗ Failed to fetch {acc} (not found in UniProt)")
                failed_fetches += 1
        except Exception as e:
            print(f"✗ Error fetching {acc}: {e}")
            failed_fetches += 1

    print(f"\nFetch summary:")
    print(f"  Successful: {successful_fetches}")
    print(f"  Failed: {failed_fetches}")
    
    if successful_fetches == 0 and len(missing_accessions) > 0:
        print("ERROR: Could not fetch any metadata from UniProt!")
        print("Possible causes:")
        print("  - Network connectivity issues")
        print("  - Invalid accession format")
        print("  - All accessions don't exist in UniProt")
        print("Run debug_metadata.py to diagnose the issue.")
    
    # Save updated metadata
    try:
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as out:
            json.dump(results, out, indent=2)
        print(f"✓ Updated metadata saved to {out_json}")
    except Exception as e:
        print(f"ERROR saving metadata: {e}")

    return results





def generate_additional_itol_files(metadata_json_path, itol_folder):
    """
    Build UniProt-based iTOL datasets to mirror CUPP group styling,
    without TREE_COLORS clades/branches, and with proper 'arrow' connections.

    What this version writes into <itol_folder>/additional:
      - Protein_description_text.txt (DATASET_TEXT, per-string color)
      - protein_domains_and_signals.txt (DATASET_DOMAINS, shows signal peptides and Pfam domains)
      - protein_supfam_domains.txt (DATASET_DOMAINS, shows SUPFAM domains)
      - GO:
          * go_terms_text_multi.txt (multi-label DATASET_TEXT, per-term colors)
          * uniprot_GO_*   color strip, text, text_single, symbols
          * uniprot_GO_arrows.txt (DATASET_CONNECTION with arrows)
      - Pfam:
          * uniprot_PFAM_* color strip, text, text_single, symbols
          * uniprot_PFAM_arrows.txt (DATASET_CONNECTION with arrows)
      - SUPFAM:
          * uniprot_SUPFAM_* color strip, text, text_single, symbols
          * uniprot_SUPFAM_arrows.txt (DATASET_CONNECTION with arrows)
      - DOI:
          * uniprot_DOI_* color strip, text, text_single, symbols
          * uniprot_DOI_arrows.txt (DATASET_CONNECTION with arrows)
      - EC Numbers:
          * uniprot_EC_* color strip, text, text_single, symbols
          * uniprot_EC_arrows.txt (DATASET_CONNECTION with arrows)
      - Annotation Scores:
          * uniprot_SCORE_* color strip, text, text_single, symbols
          * uniprot_SCORE_arrows.txt (DATASET_CONNECTION with arrows)
      - Taxonomy per-rank (species, genus, family, order, class, phylum, kingdom, superkingdom):
          * taxonomy_<rank>_* color strip, text, text_single, symbols
          * taxonomy_<rank>_arrows.txt (DATASET_CONNECTION with arrows)
      - popup_info.txt (POPUP_INFO; rich HTML popups with UniProt and AlphaFold links)

    Notes
      - We only emit IDs present in the current tree. If *_pruned.tree exists in itol_folder,
        we use its leaves; otherwise we parse the *_CUPP_gr_label.txt written earlier.
      - Arrows follow iTOL's DATASET_CONNECTION header, with DRAW_ARROWS=1 and ALIGN_TO_LABELS=1.
        Each unique string picks one anchor and draws anchor -> others in that group.  # see official template  [1](https://cran.r-project.org/web/packages/itol.toolkit/itol.toolkit.pdf)
      - POPUP_INFO uses TAB separator to safely include HTML without escaping commas.  # official template  [2](https://itol.embl.de/help/popup_info_template.txt)
    """
    import os, json, hashlib, glob

    # stable hex color per label
    def stable_color(label: str) -> str:
        import hashlib as _h
        return f"#{_h.md5(label.encode('utf-8')).hexdigest()[:6].upper()}"

    def ensure_dir(p: str):
        os.makedirs(p, exist_ok=True)
        return p

    # discover current leaf set and (optionally) order
    def discover_allowed_nodes_and_order(folder: str):
        ordered = []
        allowed = set()
        pruned = sorted(glob.glob(os.path.join(folder, "*_pruned.tree")))
        if pruned:
            try:
                from Bio import Phylo
                tree = Phylo.read(pruned[0], "newick")
                try:
                    tree.ladderize()
                except Exception:
                    pass
                ordered = [t.name for t in tree.get_terminals() if t.name]
                allowed = set(ordered)
                return allowed, ordered
            except Exception:
                pass
        # fallback: parse CUPP color strip for present leaf IDs
        cs = sorted(glob.glob(os.path.join(folder, "*_CUPP_gr_label.txt")))
        ids = []
        for path in cs:
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    in_data = False
                    for ln in fh:
                        ln = ln.rstrip("\n")
                        if not in_data:
                            if ln.strip() == "DATA":
                                in_data = True
                            continue
                        if not ln or ln.startswith("#"):
                            continue
                        parts = ln.split("\t")
                        if parts:
                            ids.append(parts[0])
            except Exception:
                continue
        return set(ids), ids

    allowed_nodes, ordered_leaves = discover_allowed_nodes_and_order(itol_folder)

    # choose a representative label per node
    def choose_representative(node2labels: Dict[str, List[str]]) -> Dict[str, str]:
        freq = {}
        for labs in node2labels.values():
            for l in set(labs):
                freq[l] = freq.get(l, 0) + 1
        rep = {}
        for node, labs in node2labels.items():
            uniq = sorted(set(labs))
            uniq.sort(key=lambda x: (-freq.get(x, 0), x))
            if uniq:
                rep[node] = uniq[0]
        return rep

    def concatenate_all_labels(node2labels: Dict[str, List[str]]) -> Dict[str, str]:
        """Concatenate all labels for each node with + separator"""
        rep = {}
        for node, labs in node2labels.items():
            uniq = sorted(set(labs))
            if uniq:
                rep[node] = "+".join(uniq)
        return rep

    # iTOL DATASET_TEXT with COMMA separator and per-item color
    # line: ID,TEXT,-1,#RRGGBB,bold,12,0
    def write_text_multilabel(path, dataset_label, node2labels, label2color):
        with open(path, "w", encoding="utf-8") as out:
            out.write("DATASET_TEXT\nSEPARATOR COMMA\n")
            out.write(f"DATASET_LABEL,{dataset_label}\nCOLOR,#000000\nDATA")
            for node, labels in sorted(node2labels.items()):
                for lab in sorted(set(labels)):
                    col = label2color.get(lab, "#000000")
                    safe = str(lab).replace(",", " ").replace("\n", " ")
                    out.write(f"\n{node},{safe},-1,{col},bold,1,0")

    # CUPP-style simple writers (you already have these helpers defined globally)
    def write_cupp_colorstrip(path, label, node2group, colors):
        write_colorstrip(path, label, node2group, colors)

    def write_cupp_text(path, label, node2group, colors):
        write_text(path, label, node2group, colors)

    def write_cupp_text_single(path, label, node2group, colors):
        write_text_single(path, label, node2group, colors)

    def write_cupp_symbols(path, label, node2group, colors):
        write_symbols(path, label, node2group, colors)

    # Proper arrow dataset (DATASET_CONNECTION). One anchor per group connects to all others.
    # Official fields and options: DRAW_ARROWS, ARROW_SIZE, ALIGN_TO_LABELS, CURVE_ANGLE, MAXIMUM_LINE_WIDTH.  [1](https://cran.r-project.org/web/packages/itol.toolkit/itol.toolkit.pdf)
    def write_connections_arrows(path, dataset_label, node2group, group_colors):
        groups = {}
        for node, g in node2group.items():
            groups.setdefault(g, []).append(node)
        anchor = {g: sorted(v)[0] for g, v in groups.items() if v}
        with open(path, "w", encoding="utf-8") as out:
            out.write("DATASET_CONNECTION\n")
            out.write("SEPARATOR COMMA\n")
            out.write(f"DATASET_LABEL,{dataset_label}\n")
            out.write("COLOR,#000000\n")
            out.write("DRAW_ARROWS,1\n")
            out.write("ARROW_SIZE,20\n")
            out.write("ALIGN_TO_LABELS,1\n")
            out.write("CURVE_ANGLE,0\n")
            out.write("MAXIMUM_LINE_WIDTH,10\n")
            out.write("DATA")
            for g, members in groups.items():
                if not members:
                    continue
                src = anchor[g]
                color = group_colors[g]
                for node in sorted(members):
                    if node == src:
                        continue
                    out.write(f"\n{src},{node},6,{color},normal,{g}")

    # POPUP_INFO writer. Use TAB to avoid comma escaping in HTML. Official template: NODE_ID, POPUP_TITLE, POPUP_CONTENT.  [2](https://itol.embl.de/help/popup_info_template.txt)
    def write_popup_info(path, records):
        with open(path, "w", encoding="utf-8") as out:
            out.write("POPUP_INFO\nSEPARATOR TAB\nDATA\n")
            for node, (title, html) in sorted(records.items()):
                # replace tabs and newlines
                t = str(title).replace("\t", " ").replace("\n", " ")
                h = str(html).replace("\t", " ").replace("\n", " ")
                out.write(f"{node}\t{t}\t{h}\n")

    # Load metadata
    print(f"Loading metadata from: {metadata_json_path}")
    if not os.path.exists(metadata_json_path):
        print(f"ERROR: Metadata file not found: {metadata_json_path}")
        return
    
    try:
        with open(metadata_json_path, "r", encoding="utf-8") as fh:
            raw_meta = json.load(fh)
        print(f"✓ Metadata loaded: {len(raw_meta)} entries")
    except Exception as e:
        print(f"ERROR loading metadata: {e}")
        return

    out_dir = ensure_dir(os.path.join(itol_folder, "additional"))
    print(f"Output directory: {out_dir}")
    
    print(f"Allowed nodes: {len(allowed_nodes)} nodes")
    if len(allowed_nodes) == 0:
        print("WARNING: No allowed nodes found! Check tree and label files.")
        print(f"Checking in folder: {itol_folder}")
        return

    # Collectors
    protein_desc = {}            # node -> description
    go_node2labels = {}          # node -> [go terms]
    pfam_node2labels = {}        # node -> [pfam ids]
    supfam_node2labels = {}      # node -> [supfam ids] 
    doi_node2labels = {}         # node -> [doi]
    ec_node2labels = {}          # node -> [ec numbers]
    annotation_scores = {}       # node -> annotation score
    tax_ranks = {rk: {} for rk in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]}
    popup_records = {}           # node -> (title, html)

    # Helper parsing of UniProt entry
    def pick_node_id(meta_key: str) -> Union[str, None]:
        if meta_key in allowed_nodes:
            return meta_key
        acc_norm = meta_key.split("_", 1)[0]
        if acc_norm in allowed_nodes:
            return acc_norm
        return None

    def get(entry, path, default=None):
        cur = entry
        try:
            for p in path:
                cur = cur[p]
            return cur
        except Exception:
            return default

    def get_protein_name(entry: dict) -> Union[str, None]:
        return get(entry, ["proteinDescription", "recommendedName", "fullName", "value"])

    def extract_go_terms(entry: dict) -> List[str]:
        out = []
        for xr in entry.get("uniProtKBCrossReferences", []):
            if xr.get("database") == "GO":
                for p in xr.get("properties", []):
                    if p.get("key") == "GoTerm" and p.get("value"):
                        out.append(p["value"])
        return out

    def extract_pfam_ids(meta_entry: dict) -> List[str]:
        lst = []
        if isinstance(meta_entry.get("pfam"), list):
            lst.extend([x for x in meta_entry["pfam"] if x])
        raw = meta_entry.get("raw") or {}
        for xr in raw.get("uniProtKBCrossReferences", []):
            if xr.get("database") == "Pfam" and xr.get("id"):
                lst.append(xr["id"])
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_supfam_ids(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        for xr in raw.get("uniProtKBCrossReferences", []):
            if xr.get("database") == "SUPFAM" and xr.get("id"):
                lst.append(xr["id"])
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_doi(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        references = raw.get("references", [])
        for ref in references:
            citation = ref.get("citation", {})
            ref_cross_refs = citation.get("citationCrossReferences", [])
            for xr in ref_cross_refs:
                if xr.get("database") == "DOI" and xr.get("id"):
                    lst.append(xr["id"])
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_ec_numbers(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        protein_desc = raw.get("proteinDescription") or {}
        
        # Check recommendedName
        recommended_name = protein_desc.get("recommendedName") or {}
        ec_numbers = recommended_name.get("ecNumbers") or []
        for ec_entry in ec_numbers:
            if ec_entry.get("value"):
                lst.append(ec_entry["value"])
        
        # Check comments for CATALYTIC ACTIVITY
        comments = raw.get("comments") or []
        for comment in comments:
            if comment.get("commentType") == "CATALYTIC ACTIVITY":
                reaction = comment.get("reaction") or {}
                ec_number = reaction.get("ecNumber")
                if ec_number:
                    lst.append(ec_number)
        
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_annotation_score(meta_entry: dict) -> int:
        raw = meta_entry.get("raw") or {}
        score = raw.get("annotationScore", 0.0)
        return int(score) if score > 0 else 0

    def extract_supfam_ids(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        for xr in raw.get("uniProtKBCrossReferences", []):
            if xr.get("database") == "SUPFAM" and xr.get("id"):
                lst.append(xr["id"])
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_doi(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        references = raw.get("references", [])
        for ref in references:
            citation = ref.get("citation", {})
            ref_cross_refs = citation.get("citationCrossReferences", [])
            for xr in ref_cross_refs:
                if xr.get("database") == "DOI" and xr.get("id"):
                    lst.append(xr["id"])
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_ec_numbers(meta_entry: dict) -> List[str]:
        lst = []
        raw = meta_entry.get("raw") or {}
        protein_desc = raw.get("proteinDescription") or {}
        
        # Check recommendedName
        recommended_name = protein_desc.get("recommendedName") or {}
        ec_numbers = recommended_name.get("ecNumbers") or []
        for ec_entry in ec_numbers:
            if ec_entry.get("value"):
                lst.append(ec_entry["value"])
        
        # Check comments for CATALYTIC ACTIVITY
        comments = raw.get("comments") or []
        for comment in comments:
            if comment.get("commentType") == "CATALYTIC ACTIVITY":
                reaction = comment.get("reaction") or {}
                ec_number = reaction.get("ecNumber")
                if ec_number:
                    lst.append(ec_number)
        
        # unique, preserve order
        seen, out = set(), []
        for x in lst:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    def extract_tax_ranks(entry: dict) -> dict:
        res = {k: None for k in tax_ranks.keys()}
        org = entry.get("organism") or {}
        species = org.get("scientificName")
        if species:
            res["species"] = species
            res["genus"] = species.split(" ", 1)[0]
        lineage = org.get("lineage") or []
        if lineage:
            res["superkingdom"] = lineage[0] if len(lineage) >= 1 else None
            genus_idx = lineage.index(res["genus"]) if res["genus"] in lineage else len(lineage) - 1
            idx = {"family": genus_idx - 1, "order": genus_idx - 2, "class": genus_idx - 3,
                   "phylum": genus_idx - 4, "kingdom": genus_idx - 5}
            for rk, i in idx.items():
                res[rk] = lineage[i] if i >= 0 else None
        return res

    def has_alphafold(entry: dict) -> bool:
        for xr in entry.get("uniProtKBCrossReferences", []):
            if xr.get("database") in {"AlphaFoldDB", "AlphaFold"}:
                return True
        return False

    # Build collectors and popups (trim overly long parts)
    processed_count = 0
    skipped_count = 0
    
    print("Processing metadata entries...")
    for meta_key, meta_entry in raw_meta.items():
        node = pick_node_id(meta_key)
        if not node:
            skipped_count += 1
            continue
            
        processed_count += 1
        if processed_count <= 5:  # Show first 5 for debugging
            print(f"  Processing: {meta_key} -> node: {node}")

        raw_entry = meta_entry.get("raw") or {}
        accession = meta_key.split("_", 1)[0]

        name = get_protein_name(raw_entry)
        if name:
            protein_desc[node] = name

        terms = extract_go_terms(raw_entry)
        if terms:
            go_node2labels.setdefault(node, []).extend(terms)

        pfams = extract_pfam_ids(meta_entry)
        if pfams:
            pfam_node2labels.setdefault(node, []).extend(pfams)

        supfams = extract_supfam_ids(meta_entry)
        if supfams:
            supfam_node2labels.setdefault(node, []).extend(supfams)

        dois = extract_doi(meta_entry)
        if dois:
            doi_node2labels.setdefault(node, []).extend(dois)

        ec_numbers = extract_ec_numbers(meta_entry)
        if ec_numbers:
            ec_node2labels.setdefault(node, []).extend(ec_numbers)

        annotation_score = extract_annotation_score(meta_entry)
        if annotation_score > 0:
            annotation_scores[node] = annotation_score

        annotation_score = extract_annotation_score(meta_entry)
        if annotation_score > 0:
            annotation_scores[node] = annotation_score

        ranks = extract_tax_ranks(raw_entry)
        for rk, val in ranks.items():
            if val:
                tax_ranks[rk][node] = val

        # popup info
        len_aa = get(raw_entry, ["sequence", "length"])
        ec = None
        ec_list = get(raw_entry, ["uniProtKBCrossReferences"]) or []
        # EC number might be in 'ecNumbers' or in catalytic activity text. Keep simple.
        if isinstance(get(raw_entry, ["ecNumbers"]), list) and get(raw_entry, ["ecNumbers"]):
            ec = "; ".join(get(raw_entry, ["ecNumbers"]))
        pe = get(raw_entry, ["proteinExistence"])
        orgn = get(raw_entry, ["organism", "scientificName"])
        uni_url = f"https://www.uniprot.org/uniprotkb/{accession}"
        af_url = f"https://alphafold.ebi.ac.uk/entry/{accession}" if has_alphafold(raw_entry) else None

        # limit long lists
        pf_short = pfams[:8] if pfams else []
        go_short = terms[:5] if terms else []

        # assemble HTML compactly; avoid very long strings
        title = f"{accession} · {name}"[:140] if name else f"{accession}"
        parts = []
        if name:
            parts.append(f"<h3 style='margin:0'>{name}</h3>")
        row1 = []
        row1.append(f"<b>Accession</b> {accession}")
        if orgn:
            row1.append(f"<b>Organism</b> {orgn}")
        if len_aa:
            row1.append(f"<b>Length</b> {len_aa} aa")
        if ec:
            row1.append(f"<b>EC</b> {ec}")
        if pe:
            row1.append(f"<b>PE</b> {pe}")
        parts.append("<p style='margin:0.2em 0'>" + " | ".join(row1) + "</p>")

        if pf_short:
            pf_links = " ".join([f"<a target='_blank' href='https://www.ebi.ac.uk/interpro/entry/pfam/{p}'>{p}</a>" for p in pf_short])
            parts.append(f"<p style='margin:0.2em 0'><b>Pfam</b> {pf_links}</p>")
        if go_short:
            safe_go = "; ".join(go_short).replace("<", "&lt;")
            parts.append(f"<p style='margin:0.2em 0'><b>GO</b> {safe_go}</p>")

        link_bits = [f"<a target='_blank' href='{uni_url}'>UniProt</a>"]
        if af_url:
            link_bits.append(f"<a target='_blank' href='{af_url}'>AlphaFold</a>")
        parts.append("<p style='margin:0.2em 0'>" + " · ".join(link_bits) + "</p>")

        html = "".join(parts)
        if len(html) > 1500:  # keep popups snappy
            html = html[:1500] + "…"
        popup_records[node] = (title, html)

    print(f"\nProcessing summary:")
    print(f"  Total metadata entries: {len(raw_meta)}")
    print(f"  Processed (matching tree): {processed_count}")
    print(f"  Skipped (not in tree): {skipped_count}")
    print(f"  Proteins with descriptions: {len(protein_desc)}")
    print(f"  Proteins with GO terms: {len(go_node2labels)}")
    print(f"  Proteins with Pfam domains: {len(pfam_node2labels)}")
    print(f"  Proteins with taxonomy: {sum(len(tax_ranks[rk]) for rk in tax_ranks)}")
    
    if processed_count == 0:
        print("WARNING: No metadata entries matched any tree nodes!")
        print("This means the accessions in metadata don't match tree leaf names.")
        print("Additional files will be empty.")

    # 1) Protein description with per-string colors (COMMA)
    if protein_desc:
        desc_path = os.path.join(out_dir, "protein_description_text.txt")
        with open(desc_path, "w", encoding="utf-8") as out:
            out.write("DATASET_TEXT\nSEPARATOR COMMA\n")
            out.write("DATASET_LABEL,Protein_description\nCOLOR,#000000\nDATA")
            for node, txt in sorted(protein_desc.items()):
                col = stable_color("DESC|" + txt)
                safe = str(txt).replace(",", " ").replace("\n", " ")
                out.write(f"\n{node},{safe},-1,{col},bold,1,0")  # official TEXT format  l TEXT + CUPP-style + ARROWS
    if go_node2labels:
        go_terms = sorted({t for ts in go_node2labels.values() for t in ts})
        go_colors = {t: stable_color("GO|" + t) for t in go_terms}
        write_text_multilabel(
            os.path.join(out_dir, "go_terms_text_multi.txt"),
            "GO_terms_text_multi",
            go_node2labels,
            go_colors
        )
        go_rep = choose_representative(go_node2labels)
        if go_rep:
            go_group_colors = {g: go_colors.get(g, stable_color("GO|" + g)) for g in sorted(set(go_rep.values()))}
            write_cupp_colorstrip(os.path.join(out_dir, "uniprot_GO_uniprot_GO_label.txt"),
                                  "uniprot_GO_label", go_rep, go_group_colors)
            write_cupp_text(os.path.join(out_dir, "uniprot_GO_uniprot_GO_text.txt"),
                            "uniprot_GO_text", go_rep, go_group_colors)
            write_cupp_text_single(os.path.join(out_dir, "uniprot_GO_uniprot_GO_text_single.txt"),
                                   "uniprot_GO_text_single", go_rep, go_group_colors)
            write_cupp_symbols(os.path.join(out_dir, "uniprot_GO_uniprot_GO_symbols.txt"),
                               "uniprot_GO_symbols", go_rep, go_group_colors)
            # arrows
            write_connections_arrows(os.path.join(out_dir, "uniprot_GO_arrows.txt"),
                                     "uniprot_GO_arrows", go_rep, go_group_colors)

    # 3) Pfam: CUPP-style + ARROWS (no multi-label text file by your request)
    if pfam_node2labels:
        pf_terms = sorted({p for ps in pfam_node2labels.values() for p in ps})
        pf_colors = {p: stable_color("PF|" + p) for p in pf_terms}
        pf_rep = concatenate_all_labels(pfam_node2labels)
        if pf_rep:
            pf_group_colors = {g: pf_colors.get(g, stable_color("PF|" + g)) for g in sorted(set(pf_rep.values()))}
            write_cupp_colorstrip(os.path.join(out_dir, "uniprot_PFAM_uniprot_PFAM_label.txt"),
                                  "uniprot_PFAM_label", pf_rep, pf_group_colors)
            write_cupp_text(os.path.join(out_dir, "uniprot_PFAM_uniprot_PFAM_text.txt"),
                            "uniprot_PFAM_text", pf_rep, pf_group_colors)
            write_cupp_text_single(os.path.join(out_dir, "uniprot_PFAM_uniprot_PFAM_text_single.txt"),
                                   "uniprot_PFAM_text_single", pf_rep, pf_group_colors)
            write_cupp_symbols(os.path.join(out_dir, "uniprot_PFAM_uniprot_PFAM_symbols.txt"),
                               "uniprot_PFAM_symbols", pf_rep, pf_group_colors)
            # arrows
            write_connections_arrows(os.path.join(out_dir, "uniprot_PFAM_arrows.txt"),
                                     "uniprot_PFAM_arrows", pf_rep, pf_group_colors)

    # 4) Taxonomy per rank: CUPP-style + ARROWS
    for rk, node2val in tax_ranks.items():
        if not node2val:
            continue
        vals = sorted(set(node2val.values()))
        colors = {v: stable_color(f"TAX|{rk}|{v}") for v in vals}
        write_cupp_colorstrip(os.path.join(out_dir, f"taxonomy_{rk}_taxonomy_{rk}_label.txt"),
                              f"taxonomy_{rk}_label", node2val, colors)
        write_cupp_text(os.path.join(out_dir, f"taxonomy_{rk}_taxonomy_{rk}_text.txt"),
                        f"taxonomy_{rk}_text", node2val, colors)
        write_cupp_text_single(os.path.join(out_dir, f"taxonomy_{rk}_taxonomy_{rk}_text_single.txt"),
                               f"taxonomy_{rk}_text_single", node2val, colors)
        write_cupp_symbols(os.path.join(out_dir, f"taxonomy_{rk}_taxonomy_{rk}_symbols.txt"),
                           f"taxonomy_{rk}_symbols", node2val, colors)
        # arrows
        write_connections_arrows(os.path.join(out_dir, f"taxonomy_{rk}_arrows.txt"),
                                 f"taxonomy_{rk}_arrows", node2val, colors)

    # 4.1) SUPFAM: CUPP-style + ARROWS
    if supfam_node2labels:
        sf_terms = sorted({p for ps in supfam_node2labels.values() for p in ps})
        sf_colors = {p: stable_color("SF|" + p) for p in sf_terms}
        sf_rep = concatenate_all_labels(supfam_node2labels)
        if sf_rep:
            sf_group_colors = {g: sf_colors.get(g, stable_color("SF|" + g)) for g in sorted(set(sf_rep.values()))}
            write_cupp_colorstrip(os.path.join(out_dir, "uniprot_SUPFAM_uniprot_SUPFAM_label.txt"),
                                  "uniprot_SUPFAM_label", sf_rep, sf_group_colors)
            write_cupp_text(os.path.join(out_dir, "uniprot_SUPFAM_uniprot_SUPFAM_text.txt"),
                            "uniprot_SUPFAM_text", sf_rep, sf_group_colors)
            write_cupp_text_single(os.path.join(out_dir, "uniprot_SUPFAM_uniprot_SUPFAM_text_single.txt"),
                                   "uniprot_SUPFAM_text_single", sf_rep, sf_group_colors)
            write_cupp_symbols(os.path.join(out_dir, "uniprot_SUPFAM_uniprot_SUPFAM_symbols.txt"),
                               "uniprot_SUPFAM_symbols", sf_rep, sf_group_colors)
            # arrows
            write_connections_arrows(os.path.join(out_dir, "uniprot_SUPFAM_arrows.txt"),
                                     "uniprot_SUPFAM_arrows", sf_rep, sf_group_colors)

    # 4.2) DOI: CUPP-style + ARROWS
    if doi_node2labels:
        doi_terms = sorted({p for ps in doi_node2labels.values() for p in ps})
        doi_colors = {p: stable_color("DOI|" + p) for p in doi_terms}
        doi_rep = concatenate_all_labels(doi_node2labels)
        if doi_rep:
            doi_group_colors = {g: doi_colors.get(g, stable_color("DOI|" + g)) for g in sorted(set(doi_rep.values()))}
            write_cupp_colorstrip(os.path.join(out_dir, "uniprot_DOI_uniprot_DOI_label.txt"),
                                  "uniprot_DOI_label", doi_rep, doi_group_colors)
            write_cupp_text(os.path.join(out_dir, "uniprot_DOI_uniprot_DOI_text.txt"),
                            "uniprot_DOI_text", doi_rep, doi_group_colors)
            write_cupp_text_single(os.path.join(out_dir, "uniprot_DOI_uniprot_DOI_text_single.txt"),
                                   "uniprot_DOI_text_single", doi_rep, doi_group_colors)
            write_cupp_symbols(os.path.join(out_dir, "uniprot_DOI_uniprot_DOI_symbols.txt"),
                               "uniprot_DOI_symbols", doi_rep, doi_group_colors)
            # arrows
            write_connections_arrows(os.path.join(out_dir, "uniprot_DOI_arrows.txt"),
                                     "uniprot_DOI_arrows", doi_rep, doi_group_colors)

    # 4.3) EC Numbers: CUPP-style + ARROWS
    if ec_node2labels:
        ec_terms = sorted({p for ps in ec_node2labels.values() for p in ps})
        ec_colors = {p: stable_color("EC|" + p) for p in ec_terms}
        ec_rep = concatenate_all_labels(ec_node2labels)
        if ec_rep:
            ec_group_colors = {g: ec_colors.get(g, stable_color("EC|" + g)) for g in sorted(set(ec_rep.values()))}
            write_cupp_colorstrip(os.path.join(out_dir, "uniprot_EC_uniprot_EC_label.txt"),
                                  "uniprot_EC_label", ec_rep, ec_group_colors)
            write_cupp_text(os.path.join(out_dir, "uniprot_EC_uniprot_EC_text.txt"),
                            "uniprot_EC_text", ec_rep, ec_group_colors)
            write_cupp_text_single(os.path.join(out_dir, "uniprot_EC_uniprot_EC_text_single.txt"),
                                   "uniprot_EC_text_single", ec_rep, ec_group_colors)
            write_cupp_symbols(os.path.join(out_dir, "uniprot_EC_uniprot_EC_symbols.txt"),
                               "uniprot_EC_symbols", ec_rep, ec_group_colors)
            # arrows
            write_connections_arrows(os.path.join(out_dir, "uniprot_EC_arrows.txt"),
                                     "uniprot_EC_arrows", ec_rep, ec_group_colors)

    # 4.4) Annotation Scores: CUPP-style + ARROWS
    if annotation_scores:
        # Use discrete annotation score values (1, 2, 3, 4, 5)
        score_terms = sorted(set(annotation_scores.values()))
        score_colors = {s: stable_color(f"SCORE|{s}") for s in score_terms}
        
        write_cupp_colorstrip(os.path.join(out_dir, "uniprot_SCORE_uniprot_SCORE_label.txt"),
                              "uniprot_SCORE_label", annotation_scores, score_colors)
        write_cupp_text(os.path.join(out_dir, "uniprot_SCORE_uniprot_SCORE_text.txt"),
                        "uniprot_SCORE_text", annotation_scores, score_colors)
        write_cupp_text_single(os.path.join(out_dir, "uniprot_SCORE_uniprot_SCORE_text_single.txt"),
                               "uniprot_SCORE_text_single", annotation_scores, score_colors)
        write_cupp_symbols(os.path.join(out_dir, "uniprot_SCORE_uniprot_SCORE_symbols.txt"),
                           "uniprot_SCORE_symbols", annotation_scores, score_colors)
        # arrows
        write_connections_arrows(os.path.join(out_dir, "uniprot_SCORE_arrows.txt"),
                                 "uniprot_SCORE_arrows", annotation_scores, score_colors)

    # 5) Domain and Signal Peptide annotation (DATASET_DOMAINS)
    def extract_domains_and_signals(entry: dict) -> List[Tuple]:
        """Extract domains and signal peptides from UniProt features and Pfam cross-references.
        Returns list of (start, end, description, feature_type, pfam_id)"""
        domains = []
        raw_entry = entry.get("raw") or {}
        features = raw_entry.get("features") or []
        
        # Get all Pfam IDs and their descriptions from cross-references
        all_pfam_ids = extract_pfam_ids(entry)
        xrefs = raw_entry.get("uniProtKBCrossReferences", [])
        pfam_descriptions = {}
        for xref in xrefs:
            if xref.get("database") == "Pfam":
                pfam_id = xref.get("id")
                entry_name = ""
                properties = xref.get("properties", [])
                for prop in properties:
                    if prop.get("key") == "EntryName":
                        entry_name = prop.get("value", "")
                        break
                if pfam_id:
                    pfam_descriptions[pfam_id] = entry_name
        
        # First, extract from features (these have positional information)
        features_pfam_ids = set()
        for feature in features:
            if feature.get("type") in ["Domain", "Signal"]:
                location = feature.get("location", {})
                start_info = location.get("start", {})
                end_info = location.get("end", {})
                
                start = start_info.get("value")
                end = end_info.get("value")
                description = feature.get("description", "")
                feature_type = feature.get("type")
                
                # Extract Pfam ID from evidences if it's a domain
                pfam_id = None
                if feature_type == "Domain":
                    evidences = feature.get("evidences", [])
                    for evidence in evidences:
                        if evidence.get("source") == "Pfam":
                            pfam_id = evidence.get("id", "")
                            features_pfam_ids.add(pfam_id)
                            break
                    
                    # If no Pfam evidence found, try to match by description
                    if not pfam_id:
                        # Try to find a Pfam domain that matches this description
                        description_lower = description.lower()
                        for pid, pdesc in pfam_descriptions.items():
                            if pid not in features_pfam_ids:  # Not already assigned
                                # Check if descriptions match (e.g., "OmpA-like" matches "OmpA")
                                if (description_lower in pdesc.lower() or 
                                    pdesc.lower() in description_lower or
                                    (description_lower == "ompa-like" and pdesc.lower() == "ompa")):
                                    pfam_id = pid
                                    features_pfam_ids.add(pfam_id)
                                    break
                
                if start is not None and end is not None:
                    domains.append((start, end, description, feature_type, pfam_id or ""))
        
        # Add any remaining missing Pfam domains without positional info
        missing_pfam_ids = set(all_pfam_ids) - features_pfam_ids  
        for pfam_id in missing_pfam_ids:
            if pfam_id:  # Skip empty IDs
                entry_name = pfam_descriptions.get(pfam_id, "")
                description = entry_name if entry_name else pfam_id
                # Add as domain with no position (start=0, end=0 indicates no position)
                domains.append((0, 0, description, "Domain", pfam_id))
        
        return domains

    def extract_supfam_domains(entry: dict) -> List[Tuple]:
        """Extract SUPFAM domains from UniProt cross-references.
        Returns list of (start, end, description, feature_type, supfam_id)"""
        domains = []
        raw_entry = entry.get("raw") or {}
        
        # Get SUPFAM domains from cross-references
        xrefs = raw_entry.get("uniProtKBCrossReferences", [])
        for xref in xrefs:
            if xref.get("database") == "SUPFAM":
                supfam_id = xref.get("id")
                entry_name = ""
                properties = xref.get("properties", [])
                for prop in properties:
                    if prop.get("key") == "EntryName":
                        entry_name = prop.get("value", "")
                        break
                
                if supfam_id:
                    description = entry_name if entry_name else supfam_id
                    # SUPFAM domains typically don't have position info in UniProt, so use 0,0
                    domains.append((0, 0, description, "SUPFAM", supfam_id))
        
        return domains

    # DATASET_DOMAINS writer
    def write_domains_file(path, dataset_label, node_domains, protein_lengths):
        """Write iTOL DATASET_DOMAINS file with signal peptides and Pfam domains"""
        with open(path, "w", encoding="utf-8") as out:
            out.write("DATASET_DOMAINS\n")
            out.write("SEPARATOR TAB\n")
            out.write(f"DATASET_LABEL\t{dataset_label}\n")
            out.write("COLOR\t#ff0000\n")
            out.write("DATA\n")
            
            for node, domains in sorted(node_domains.items()):
                if not domains:
                    continue
                
                protein_length = protein_lengths.get(node, 1000)  # default length if unknown
                domain_parts = []
                
                for start, end, description, feature_type, pfam_id in domains:
                    # Skip domains with no positional information
                    if start == 0 and end == 0:
                        continue
                        
                    if feature_type == "Signal":
                        # Signal peptides - green color
                        label = "Signal"
                        color = "#90EE90"  # Light green
                        domain_parts.append(f"RE|{start}|{end}|{color}|{label}")
                    elif feature_type == "Domain" and pfam_id:
                        # Pfam domains - blue color, show Pfam ID
                        color = stable_color(pfam_id)  # Use stable color per Pfam family
                        label = f"{pfam_id}"
                        if description and description != pfam_id:
                            label = f"{pfam_id} ({description})"
                        domain_parts.append(f"HH|{start}|{end}|{color}|{label}")
                    elif feature_type == "Domain":
                        # Other domains without Pfam ID - purple color  
                        color = "#9370DB"  # Medium slate blue
                        label = description or "Domain"
                        domain_parts.append(f"HH|{start}|{end}|{color}|{label}")
                
                if domain_parts:
                    domain_string = "\t".join(domain_parts)
                    out.write(f"{node}\t{protein_length}\t{domain_string}\n")

    # SUPFAM DATASET_DOMAINS writer
    def write_supfam_domains_file(path, dataset_label, node_domains, protein_lengths):
        """Write iTOL DATASET_DOMAINS file with SUPFAM domains"""
        with open(path, "w", encoding="utf-8") as out:
            out.write("DATASET_DOMAINS\n")
            out.write("SEPARATOR TAB\n")
            out.write(f"DATASET_LABEL\t{dataset_label}\n")
            out.write("COLOR\t#ff0000\n")
            out.write("DATA\n")
            
            for node, domains in sorted(node_domains.items()):
                if not domains:
                    continue
                
                protein_length = protein_lengths.get(node, 1000)  # default length if unknown
                domain_parts = []
                
                for start, end, description, feature_type, supfam_id in domains:
                    if feature_type == "SUPFAM" and supfam_id:
                        # SUPFAM domains - use specific colors for superfamilies
                        color = stable_color(supfam_id)  # Use stable color per SUPFAM family
                        label = f"{supfam_id}"
                        if description and description != supfam_id:
                            label = f"{supfam_id} ({description})"
                        # Since SUPFAM typically doesn't have positions, show as full-length rectangle
                        domain_parts.append(f"RE|1|{protein_length}|{color}|{label}")
                
                if domain_parts:
                    domain_string = "\t".join(domain_parts)
                    out.write(f"{node}\t{protein_length}\t{domain_string}\n")

    # Collect domain information 
    node_domains = {}
    protein_lengths = {}
    
    for meta_key, meta_entry in raw_meta.items():
        node = pick_node_id(meta_key)
        if not node:
            continue
        
        raw_entry = meta_entry.get("raw") or {}
        
        # Get protein length
        length = get(raw_entry, ["sequence", "length"])
        if length:
            protein_lengths[node] = length
        
        # Extract domains and signal peptides
        domains = extract_domains_and_signals(meta_entry)
        if domains:
            node_domains[node] = domains
    
    # Write domains file if we have domain data
    if node_domains:
        domains_path = os.path.join(out_dir, "protein_domains_and_signals.txt")
        write_domains_file(domains_path, "Domains_and_Signals", node_domains, protein_lengths)

    # Collect SUPFAM domain information 
    node_supfam_domains = {}
    
    for meta_key, meta_entry in raw_meta.items():
        node = pick_node_id(meta_key)
        if not node:
            continue
        
        # Extract SUPFAM domains
        supfam_domains = extract_supfam_domains(meta_entry)
        if supfam_domains:
            node_supfam_domains[node] = supfam_domains
    
    # Write SUPFAM domains file if we have domain data
    if node_supfam_domains:
        supfam_domains_path = os.path.join(out_dir, "protein_supfam_domains.txt")
        write_supfam_domains_file(supfam_domains_path, "SUPFAM_Domains", node_supfam_domains, protein_lengths)

    # 6) POPUP_INFO with rich HTML and AlphaFold links when available
    if popup_records:
        write_popup_info(os.path.join(out_dir, "popup_info.txt"), popup_records)

    print(f"### UniProt iTOL files written to: {out_dir}")
    if node_domains:
        print(f"### Domain annotation file: protein_domains_and_signals.txt")
    if node_supfam_domains:
        print(f"### SUPFAM domains file: protein_supfam_domains.txt")
    if supfam_node2labels:
        print(f"### SUPFAM data files: uniprot_SUPFAM_*")
    if doi_node2labels:
        print(f"### DOI data files: uniprot_DOI_*")
    if ec_node2labels:
        print(f"### EC number data files: uniprot_EC_*")
    if annotation_scores:
        print(f"### Annotation score data files: uniprot_SCORE_*")










def clean_acc(label: str) -> str:
    """Return accession as first token and without any :suffix decorations."""
    if label is None:
        return None
    return label.split(" ", 1)[0].split(":", 1)[0].strip()


def hex_color() -> str:
    """Generate a readable random hex color while avoiding pale yellow."""
    while True:
        r, g, b = random.randrange(256), random.randrange(256), random.randrange(256)
        if not (r > 200 and g > 200):
            return "#{:02X}{:02X}{:02X}".format(r, g, b)


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


# ------------------------------
# Load CUPP Pool
# ------------------------------

def load_pool(pool_path: str, setting: Union[str, None]):
    with open(pool_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict) or not data:
        raise SystemExit(f"Invalid or empty pool JSON: {pool_path}")

    if setting:
        if setting not in data:
            raise SystemExit(f"Setting '{setting}' not found in pool. Available: {list(data.keys())}")
        key = setting
    else:
        key = list(data.keys())[0]

    block = data[key]
    if "meta" not in block or "meta_categories" not in block:
        raise SystemExit("Pool JSON missing 'meta' or 'meta_categories'.")

    return key, block["meta"], block["meta_categories"]


# ------------------------------
# iTOL dataset writers
# ------------------------------
###############################################################################
# iTOL DATASET WRITERS (FINAL, FIXED)
###############################################################################

def write_colorstrip(path, dataset_label, node2group, group_colors):
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{dataset_label}\nCOLOR\t#000000\nDATA")
        for node, g in node2group.items():
            out.write(f"\n{node}\t{group_colors[g]}\t{g}")


def write_text(path, dataset_label, node2group, group_colors):
    """
    Regular floating TEXT labels (all leaves).
    """
    with open(path, "w") as out:
        out.write("DATASET_TEXT\nSEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{dataset_label}\nCOLOR,#000000\nDATA")
        for node, g in node2group.items():
            out.write(f"\n{node},{g},-1,{group_colors[g]},bold,1,0")


def write_text_single(path, dataset_label, node2group, group_colors):
    """
    SINGLE representative text label per group.
    Choose the FIRST occurrence of each group.
    Others get NO label.
    """
    # pick one representative leaf per group
    group_rep = {}
    for node, g in node2group.items():
        if g not in group_rep:
            group_rep[g] = node

    with open(path, "w") as out:
        out.write("DATASET_TEXT\nSEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{dataset_label}\nCOLOR,#000000\nDATA")
        for g, node in group_rep.items():
            out.write(f"\n{node},{g},-1,{group_colors[g]},bold,3,0")


def write_symbols(path, dataset_label, node2group, group_colors):
    """
    Filled circle symbols (all leaves).
    internal = 1 → filled
    """
    with open(path, "w") as out:
        out.write("DATASET_SYMBOL\nSEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{dataset_label}\nCOLOR,#000000\nDATA")
        for node, g in node2group.items():
            out.write(f"\n{node},1,8,{group_colors[g]},1,1,1")

def write_arrows(path, dataset_label, node2group, group_colors):
    """
    Create a valid iTOL DATASET_CONNECTION file with arrowheads.
    For each distinct group in node2group, choose one anchor leaf and connect
    anchor -> every other leaf in that same group.

    Output format (per iTOL template):
      Header:
        DATASET_CONNECTION
        SEPARATOR COMMA
        DATASET_LABEL,<label>
        COLOR,<fallback color>
        DRAW_ARROWS,1
        ARROW_SIZE,20
        ALIGN_TO_LABELS,1
        CURVE_ANGLE,0
        MAXIMUM_LINE_WIDTH,10
        DATA

      Data line:
        NODE1,NODE2,WIDTH,COLOR,STYLE,LABEL

    Reference: iTOL DATASET_CONNECTION template and options.  # [1](https://itol.embl.de/help/dataset_connections_template.txt)
    """
    # group -> [nodes]
    groups = {}
    for node, g in node2group.items():
        groups.setdefault(g, []).append(node)

    # choose one anchor per group: pick the lexicographically first to be deterministic
    group_anchor = {g: sorted(members)[0] for g, members in groups.items() if members}

    with open(path, "w", encoding="utf-8") as out:
        out.write("DATASET_CONNECTION\n")
        out.write("SEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{dataset_label}\n")
        out.write("COLOR,#000000\n")
        out.write("DRAW_ARROWS,1\n")         # enable arrowheads on destination end  # [1](https://itol.embl.de/help/dataset_connections_template.txt)
        out.write("ARROW_SIZE,20\n")         # arrow size in pixels                   # [1](https://itol.embl.de/help/dataset_connections_template.txt)
        out.write("ALIGN_TO_LABELS,1\n")     # start/end near leaf labels             # [1](https://itol.embl.de/help/dataset_connections_template.txt)
        out.write("CURVE_ANGLE,0\n")         # straight lines                         # [1](https://itol.embl.de/help/dataset_connections_template.txt)
        out.write("MAXIMUM_LINE_WIDTH,10\n") # global scaling for widths              # [1](https://itol.embl.de/help/dataset_connections_template.txt)
        out.write("DATA")

        for g, members in groups.items():
            if not members:
                continue
            anchor = group_anchor[g]
            color = group_colors[g]
            # connect anchor -> others in the group
            for node in sorted(members):
                if node == anchor:
                    continue
                # NODE1,NODE2,WIDTH,COLOR,STYLE,LABEL
                out.write(f"\n{anchor},{node},6,{color},normal,{g}")


def write_clades_range_by_endpoints(output_path: str,
                                    tree_path: str,
                                    node2group: dict,
                                    group_colors: dict):
    """
    Write iTOL TREE_COLORS so that each group's internal lines are colored
    using a 'range' between the first and last member of the group along the
    current leaf order of the tree.

    Format written:
      TREE_COLORS
      SEPARATOR TAB
      DATA
      <leafA>|<leafB>    range    <#hex>    <label>     # when group size >= 2
      <leaf>             branch   <#hex>    <label>     # when group size == 1
    """
    # Load tree and get current ordered list of leaf names
    tree = Phylo.read(tree_path, "newick")
    # Optional: stabilize ordering for reproducibility
    try:
        tree.ladderize()
    except Exception:
        pass

    ordered_leaves = [term.name for term in tree.get_terminals()]
    pos = {name: i for i, name in enumerate(ordered_leaves)}

    # Build group -> members that actually exist in this tree
    groups = {}
    for node, g in node2group.items():
        if node in pos:
            groups.setdefault(g, []).append(node)

    with open(output_path, "w", encoding="utf-8") as out:
        out.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")

        for g, members in sorted(groups.items()):
            if not members:
                continue
            col = group_colors[g]

            # Sort this group's members by their position in the tree
            members_sorted = sorted(members, key=lambda n: pos[n])

            if len(members_sorted) == 1:
                # Single member: color just that branch
                out.write(f"{members_sorted[0]}\tbranch\t{col}\t{g}\n")
            else:
                # Use first and last member as endpoints of the colored range
                first_leaf = members_sorted[0]
                last_leaf = members_sorted[-1]
                out.write(f"{first_leaf}|{last_leaf}\trange\t{col}\t{g}\n")

# ------------------------------
# Tree pruning with BioPython
# ------------------------------

def prune_tree(tree_path: str, keep_set_clean: Set[str], out_path: str):
    """
    Prune the tree so that only leaves whose CLEAN accession is in keep_set_clean remain.
    Branch lengths are preserved. No regex is used.
    """
    tree = Phylo.read(tree_path, "newick")

    to_remove = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            c = clean_acc(clade.name)
            if c not in keep_set_clean:
                to_remove.append(clade)

    for leaf in to_remove:
        try:
            tree.prune(leaf)
        except ValueError:
            pass

    Phylo.write(tree, out_path, "newick")


# ------------------------------
# Main
# ------------------------------

def main():
    
    p = argparse.ArgumentParser(description="Generate iTOL datasets for CUPP groups and prune the tree to valid members.")
    p.add_argument("-pool", required=True, help="Path to <FAM>_CUPPpool.json")
    p.add_argument("-tree", required=True, help="Path to the original Newick tree")
    p.add_argument("-name", default="", help="Dataset name prefix, e.g. GH29")
    p.add_argument("-outdir", default="", help="Output directory, default next to the pool file under itol/")
    p.add_argument("-setting", default="", help="Optional setting key in the pool JSON")
    p.add_argument("-fasta", required=True, help="Path to the cleaned FASTA file")
    args = p.parse_args()

    # Load pool
    setting, meta, meta_cat = load_pool(args.pool, args.setting if args.setting else None)

    # Output folder
    base = os.path.dirname(os.path.dirname(os.path.abspath(args.pool)))
    name = args.name if args.name else os.path.basename(args.pool).replace("_CUPPpool.json", "")
    outdir = ensure_dir(args.outdir if args.outdir else os.path.join(base, "itol"))
    ensure_dir(os.path.join(outdir, name))

    # Load tree and build mapping from CLEAN accession to actual leaf name
    tree = Phylo.read(args.tree, "newick")
    clean_to_actual = {}
    for clade in tree.find_clades():
        if clade.is_terminal() and clade.name:
            c = clean_acc(clade.name)
            if c:
                clean_to_actual[c] = clade.name

    if not clean_to_actual:
        raise SystemExit("Tree has no named leaves to map against.")

    tree_clean_set = set(clean_to_actual.keys())

    # Build clean accession -> group using only leaves present in the tree
    if "Accession" not in meta_cat:
        raise SystemExit("Pool meta_categories does not contain 'Accession'.")

    acc_idx = meta_cat.index("Accession")
    acc_clean_to_group = {}

    for group, fields in meta.items():
        acc_dict = fields[acc_idx]   # dict of accession -> score-like value
        for acc in acc_dict.keys():
            c = clean_acc(acc)
            if c in tree_clean_set:          # keep only those that exist in the tree
                acc_clean_to_group[c] = group

    if not acc_clean_to_group:
        raise SystemExit("No CUPP accessions overlap with tree leaves. Nothing to label.")

    # Convert to node-name -> group using the tree's actual leaf labels
    node2group = {clean_to_actual[c]: g for c, g in acc_clean_to_group.items()}

    # Assign colors per group
    groups = sorted(set(acc_clean_to_group.values()))
    group_colors = {g: hex_color() for g in groups}

    # Write datasets
    write_colorstrip(os.path.join(outdir+f"/{name}", f"{name}_CUPP_gr_label.txt"),
                     f"{name}_CUPP_gr_label", node2group, group_colors)

    write_text(os.path.join(outdir+f"/{name}", f"{name}_CUPP_gr_text.txt"),
               f"{name}_CUPP_gr_text", node2group, group_colors)

    write_text_single(os.path.join(outdir+f"/{name}", f"{name}_CUPP_gr_text_single.txt"),
                      f"{name}_CUPP_gr_text_single", node2group, group_colors)

    write_arrows(os.path.join(outdir+f"/{name}", f"{name}_CUPP_gr_arrow.txt"),
                 f"{name}_CUPP_gr_arrow", node2group, group_colors)

    write_symbols(os.path.join(outdir+f"/{name}", f"{name}_CUPP_gr_symbols.txt"),
                  f"{name}_CUPP_gr_symbols", node2group, group_colors)

    # Load existing metadata OR download if missing
    print("\n--- UniProt Metadata Processing ---")
    accessions = load_accessions_from_fasta(args.fasta)
    print(f"Extracted {len(accessions)} accessions from FASTA file")
    if not accessions:
        print("WARNING: No accessions found in FASTA file! Check header format.")
        print("Expected format: >prefix|ACCESSION|suffix")
        print("Skipping additional iTOL files generation.")
    else:
        print(f"Sample accessions: {accessions[:5]}{'...' if len(accessions) > 5 else ''}")
        
        try:
            metadata = load_or_fetch_metadata(accessions, "CUPP/metadata_output.json")
            print(f"Metadata loaded/fetched successfully: {len(metadata)} entries")
            
            print("Generating additional iTOL annotation files...")
            generate_additional_itol_files("CUPP/metadata_output.json", f"{outdir}/{name}")
            print("✓ Additional iTOL files generated successfully")
            
        except Exception as e:
            print(f"ERROR during metadata processing: {e}")
            print("Additional iTOL files may not be generated correctly.")
            import traceback
            traceback.print_exc()

    # Prune the tree first so coloring matches the final leaf set
    pruned_tree = os.path.join(outdir, f"{name}/{name}_pruned.tree")
    prune_tree(args.tree, set(acc_clean_to_group.keys()), pruned_tree)

    # Create TREE_COLORS using first and last member as endpoints per group
    write_clades_range_by_endpoints(
        os.path.join(outdir, f"{name}/{name}_CUPP_gr_clades.txt"),
        pruned_tree,
        node2group,
        group_colors
    )

    print("### CUPPvisualization completed")
    print(f"### Datasets written in: {outdir}")
    print(f"### Pruned tree: {pruned_tree}")
    print(f"### Groups: {len(groups)}")
    print(f"### Valid accessions: {len(acc_clean_to_group)}")


if __name__ == "__main__":
    main()