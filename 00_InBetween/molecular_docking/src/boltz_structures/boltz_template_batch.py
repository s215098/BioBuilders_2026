#!/usr/bin/env python3
"""
Generic batch driver for Boltz-2.

Takes a YAML *template* with a `{seq}` placeholder (and optionally `{accession}`),
substitutes per-sequence, and runs Boltz over the whole FASTA. Resume-safe.
Output structure + metrics are identical to the dedicated batch scripts:

    <out_dir>/<accession>/boltz_input.yaml
    <out_dir>/<accession>/boltz_out/...
    <out_dir>/structures_pdb/<accession>.pdb
    <out_dir>/boltz_<tag>_summary.csv

The CSV captures whichever metrics Boltz produced (confidence_*, affinity_*).
If your template has an `affinity:` property, affinity columns will be populated;
otherwise they'll be blank.

Examples:
  # Heme + NNBT, affinity on NNBT:
  python boltz_template_batch.py -f boltz_input_100.fasta -t templates/heme_nnbt.yaml -o outputs_nnbt

  # 4 copper ions, no substrate:
  python boltz_template_batch.py -f cu_enzymes.fasta -t templates/4cu.yaml -o outputs_4cu

  # Custom: protein + heme + your own substrate, affinity on substrate:
  python boltz_template_batch.py -f sequences.fasta -t my_template.yaml -o outputs_custom

Template format (free-form YAML, anywhere `{seq}` or `{accession}` appears
is substituted per sequence):

    sequences:
      - protein:
          id: A
          sequence: "{seq}"
      - ligand:
          id: B
          ccd: HEM
      - ligand:
          id: C
          smiles: 'CC1=CC=C(C=C1)N(CC(C)O)CC(C)O'
    properties:
      - affinity:
          binder: C
"""

import os
import re
import json
import csv
import glob
import shutil
import argparse
import subprocess
from pathlib import Path


def log(msg, tag="BOLTZ"):
    print(f"[{tag}] {msg}", flush=True)


def parse_fasta(path):
    """Yield (accession, sequence). Accession = first whitespace token after '>'."""
    acc, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(seq)
                acc = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if acc is not None:
        yield acc, "".join(seq)


def render_yaml(template_text, sequence, accession, out_path):
    """Replace {seq} and {accession} in the template and write to out_path.

    Uses plain str.replace (not .format) so YAML braces in the template are safe.
    """
    rendered = template_text.replace("{seq}", sequence).replace("{accession}", accession)
    with open(out_path, "w") as f:
        f.write(rendered)
    return out_path


def run_boltz(yaml_path, boltz_out):
    subprocess.run(
        [
            "boltz", "predict", yaml_path,
            "--out_dir", boltz_out,
            "--use_msa_server",
            "--output_format", "pdb",
            "--no_kernels",
            "--override",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )


def to_pdb(struct_path, dest, tag="BOLTZ"):
    """Copy struct to dest as PDB. If it's a .cif, convert it (gemmi preferred)."""
    if struct_path.lower().endswith(".pdb"):
        shutil.copyfile(struct_path, dest)
        return True
    try:
        import gemmi
        st = gemmi.read_structure(struct_path)
        st.setup_entities()
        st.write_pdb(dest)
        return True
    except Exception:
        try:
            from Bio.PDB import MMCIFParser, PDBIO
            s = MMCIFParser(QUIET=True).get_structure("x", struct_path)
            io = PDBIO()
            io.set_structure(s)
            io.save(dest)
            return True
        except Exception as e:
            log(f"  could not convert {struct_path} to PDB: {e}", tag)
            return False


def first(pattern):
    hits = glob.glob(pattern, recursive=True)
    return hits[0] if hits else None


def collect_metrics(boltz_out):
    """Return (metrics_dict, structure_path). Captures whichever keys exist."""
    metrics = {}

    conf_path = first(f"{boltz_out}/**/confidence_*.json")
    if conf_path:
        with open(conf_path) as f:
            c = json.load(f)
        for k in (
            "confidence_score", "ptm", "iptm", "ligand_iptm",
            "protein_iptm", "complex_plddt", "complex_iplddt",
        ):
            if k in c:
                metrics[k] = c[k]

    aff_path = first(f"{boltz_out}/**/affinity_*.json")
    if aff_path:
        with open(aff_path) as f:
            a = json.load(f)
        for k in ("affinity_pred_value", "affinity_probability_binary"):
            if k in a:
                metrics[k] = a[k]

    struct = first(f"{boltz_out}/**/*.cif") or first(f"{boltz_out}/**/*.pdb")
    return metrics, struct


# Always written in the same order; missing values just stay blank.
FIELDS = [
    "accession", "status", "seq_len",
    "confidence_score", "ptm", "iptm", "ligand_iptm",
    "protein_iptm", "complex_plddt", "complex_iplddt",
    "affinity_pred_value", "affinity_probability_binary",
    "pdb_path", "structure_path",
]


def validate_template(text):
    """Light sanity check: must contain {seq} and have a sequences: block."""
    if "{seq}" not in text:
        raise SystemExit(
            "Template has no `{seq}` placeholder. Put `sequence: \"{seq}\"` "
            "where the FASTA sequence should be substituted."
        )
    if not re.search(r"^\s*sequences\s*:", text, re.MULTILINE):
        raise SystemExit("Template has no top-level `sequences:` block.")


def main():
    p = argparse.ArgumentParser(
        description="Generic batch Boltz-2 driver with a YAML template",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("-f", "--fasta",    required=True, help="Input FASTA")
    p.add_argument("-t", "--template", required=True, help="YAML template with {seq} placeholder")
    p.add_argument("-o", "--out_dir",  required=True, help="Root output directory")
    p.add_argument("--csv", default=None,
                   help="Summary CSV path (default: <out_dir>/boltz_<template_stem>_summary.csv)")
    p.add_argument("--resume", action="store_true",
                   help="Skip accessions that already have a structure file")
    p.add_argument("--pdb-dir", default=None,
                   help="Folder for <accession>.pdb copies (default: <out_dir>/structures_pdb)")
    p.add_argument("--tag", default=None,
                   help="Log tag (default: from template filename, e.g. BOLTZ-HEME-NNBT)")
    args = p.parse_args()

    template_text = Path(args.template).read_text()
    validate_template(template_text)

    stem = Path(args.template).stem  # e.g. heme_nnbt
    tag = args.tag or f"BOLTZ-{stem.upper().replace('_','-')}"

    os.makedirs(args.out_dir, exist_ok=True)
    csv_path = args.csv or os.path.join(args.out_dir, f"boltz_{stem}_summary.csv")
    pdb_dir = args.pdb_dir or os.path.join(args.out_dir, "structures_pdb")
    os.makedirs(pdb_dir, exist_ok=True)

    records = list(parse_fasta(args.fasta))
    log(f"Parsed {len(records)} sequences from {args.fasta}.", tag)
    log(f"Template: {args.template}", tag)
    log(f"Output:   {args.out_dir}", tag)

    fresh = not (args.resume and os.path.exists(csv_path))
    csv_f = open(csv_path, "a", newline="")
    writer = csv.DictWriter(csv_f, fieldnames=FIELDS, extrasaction="ignore")
    if fresh:
        writer.writeheader()
        csv_f.flush()

    for i, (acc, seq) in enumerate(records, 1):
        target_dir = os.path.join(args.out_dir, acc)
        boltz_out = os.path.join(target_dir, "boltz_out")
        os.makedirs(target_dir, exist_ok=True)

        print(f"\n{'='*60}\n [{i}/{len(records)}] {acc}  (len={len(seq)})\n{'='*60}")

        if args.resume and (first(f"{boltz_out}/**/*.cif") or first(f"{boltz_out}/**/*.pdb")):
            log("Already done, collecting existing metrics (resume).", tag)
            metrics, struct = collect_metrics(boltz_out)
            pdb_dest = os.path.join(pdb_dir, f"{acc}.pdb")
            ok = to_pdb(struct, pdb_dest, tag) if struct else False
            row = {"accession": acc, "status": "ok", "seq_len": len(seq),
                   "pdb_path": pdb_dest if ok else "", "structure_path": struct or ""}
            row.update(metrics)
            writer.writerow(row)
            csv_f.flush()
            continue

        try:
            yaml_path = os.path.join(target_dir, "boltz_input.yaml")
            render_yaml(template_text, seq, acc, yaml_path)
            log("Running Boltz-2...", tag)
            run_boltz(yaml_path, boltz_out)
            metrics, struct = collect_metrics(boltz_out)
            if not metrics:
                raise RuntimeError("Boltz finished but no confidence/affinity JSON found")
            pdb_dest = os.path.join(pdb_dir, f"{acc}.pdb")
            ok = to_pdb(struct, pdb_dest, tag) if struct else False
            row = {"accession": acc, "status": "ok", "seq_len": len(seq),
                   "pdb_path": pdb_dest if ok else "", "structure_path": struct or ""}
            row.update(metrics)
            log(
                "confidence={c} ptm={p} ligand_iptm={l} affinity={a}".format(
                    c=metrics.get("confidence_score"),
                    p=metrics.get("ptm"),
                    l=metrics.get("ligand_iptm"),
                    a=metrics.get("affinity_pred_value"),
                ), tag,
            )
        except Exception as e:
            log(f"FAILED: {e}", tag)
            row = {"accession": acc, "status": f"error: {e}", "seq_len": len(seq)}

        writer.writerow(row)
        csv_f.flush()

    csv_f.close()
    log(f"Done. Summary written to {csv_path}", tag)


if __name__ == "__main__":
    main()
