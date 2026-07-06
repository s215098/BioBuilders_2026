#!/usr/bin/env python3
"""End-to-end pocket detection, filtering, docking, and pose analysis.

Stages
------
  0. prepare   — convert PDB → PDBQT for any receptor without a matching .pdbqt
                 (opt-in via --prepare; uses prepare_receptor.py)
  1. fpocket   — run fpocket on every PDB in input_dir (via run_fpocket.py)
  2. filter    — keep structures whose Pocket 1 passes loose metric thresholds
  3. docking   — targeted Vina docking on passing structures (auto-computes box
                 from pocket1_vert.pqr + heme Fe in the receptor PDB)
  4. analysis  — per-pose orientation check: ring toward heme, OHs outward

Only structures that have a matching .pdbqt in input_dir proceed to docking.
Stages can be skipped with --skip-fpocket / --skip-filter / --skip-docking.

Examples
--------
    # Full pipeline including receptor preparation, NNBT ligand
    python3 pipeline.py inputs/ -l ligands/NNBT.sdf --prepare

    # Skip fpocket (outputs/ already populated) and tighten filter
    python3 pipeline.py inputs/ -l ligands/NNBT.sdf --skip-fpocket \\
        --druggability 0.4 --volume 500

    # Only run filtering to see what would pass
    python3 pipeline.py inputs/ --skip-fpocket --skip-docking
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.insert(0, str(SCRIPT_DIR))

from pocket_box import compute_box, find_heme_fe
from filter_pockets import load_all_metrics, apply_filters, find_heme_pocket
from analyze_poses import analyze_file, print_results
from prepare_receptor import convert, prepare_directory


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, label):
    """Run a subprocess command, print live output, raise on failure."""
    print(f"\n[{label}] {' '.join(str(c) for c in cmd)}")
    result = subprocess.run([str(c) for c in cmd], text=True,
                            capture_output=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {result.returncode})")


def fpocket_output_dir(outputs_dir, stem):
    return Path(outputs_dir) / f"{stem}_out"


def pocket_pqr(outputs_dir, stem, pocket_n=1):
    return fpocket_output_dir(outputs_dir, stem) / "pockets" / f"pocket{pocket_n}_vert.pqr"


def find_receptor_pdb(stem, input_dir, outputs_dir):
    """Return a PDB file suitable for heme Fe detection.

    Prefers the original input PDB. Falls back to fpocket's own output PDB
    (<stem>_out/<stem>_out.pdb), which fpocket always writes in PDB format
    even when the original input was a .cif file.
    """
    pdb = Path(input_dir) / f"{stem}.pdb"
    if pdb.exists():
        return pdb
    fallback = Path(outputs_dir) / f"{stem}_out" / f"{stem}_out.pdb"
    if fallback.exists():
        return fallback
    raise FileNotFoundError(
        f"No PDB found for {stem} in {input_dir} or {outputs_dir}/{stem}_out/"
    )


def derive_box(stem, input_dir, outputs_dir, buffer, heme_max_dist, pocket_n=1):
    """Compute Vina box for the selected pocket, with optional heme clipping."""
    pqr = pocket_pqr(outputs_dir, stem, pocket_n)
    if not pqr.exists():
        raise FileNotFoundError(f"PQR not found: {pqr}")

    try:
        receptor_pdb = find_receptor_pdb(stem, input_dir, outputs_dir)
        heme_ref = find_heme_fe(str(receptor_pdb))
    except (FileNotFoundError, ValueError) as e:
        print(f"    [WARN] {e} — box computed without heme clipping")
        heme_ref = None

    center, size = compute_box(
        str(pqr),
        buffer=buffer,
        max_dist=heme_max_dist if heme_ref is not None else None,
        heme_ref=heme_ref,
    )
    return center.tolist(), size.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("input_dir",
                   help="Directory containing .pdb (and optionally .pdbqt) files")
    p.add_argument("-l", "--ligand", default=str(SCRIPT_DIR / "ligands" / "NNBT.sdf"),
                   help="Ligand SDF file (default: ligands/NNBT.sdf)")

    dirs = p.add_argument_group("directories")
    dirs.add_argument("--outputs-dir", default=str(SCRIPT_DIR / "outputs"),
                      help="Where fpocket _out/ folders live (default: outputs/)")
    dirs.add_argument("--docking-out", default=str(SCRIPT_DIR / "docking_out"),
                      help="Where Vina poses are saved (default: docking_out/)")

    skip = p.add_argument_group("stages")
    skip.add_argument("--prepare", action="store_true",
                      help="Stage 0: convert PDB → PDBQT for any receptor that "
                           "lacks a matching .pdbqt in input_dir (uses prepare_receptor.py)")
    skip.add_argument("--keep-water", action="store_true",
                      help="When --prepare is active, retain crystallographic waters")
    skip.add_argument("--skip-fpocket", action="store_true",
                      help="Skip running fpocket (assume outputs/ is already populated)")
    skip.add_argument("--skip-filter", action="store_true",
                      help="Skip pocket filtering (dock all structures)")
    skip.add_argument("--skip-docking", action="store_true",
                      help="Skip docking (only run fpocket + filter)")

    box = p.add_argument_group("docking box")
    box.add_argument("--heme-max-dist", type=float, default=12.0,
                     help="Clip box to this Å from heme Fe (default: 12)")
    box.add_argument("--buffer", type=float, default=6.0,
                     help="Box buffer per axis in Å (default: 6.0)")
    box.add_argument("--box-size", type=float, nargs=3, default=None,
                     metavar=("X", "Y", "Z"),
                     help="Override automatic box size with fixed dimensions")

    vina = p.add_argument_group("vina")
    vina.add_argument("--exhaustiveness", type=int, default=32,
                      help="Vina exhaustiveness (default: 32)")
    vina.add_argument("--n-poses", type=int, default=10,
                      help="Poses per docking run (default: 10)")

    filt = p.add_argument_group("filter thresholds (heme-proximal pocket)")
    filt.add_argument("--fe-proximity", type=float, default=10.0,
                      help="Max distance (Å) from any alpha sphere to heme Fe "
                           "for a pocket to be considered heme-proximal (default: 10.0)")
    filt.add_argument("--druggability", type=float, default=0.2,
                      help="Min druggability score (default: 0.2)")
    filt.add_argument("--volume", type=float, default=500.0,
                      help="Min volume Å³ (default: 500)")
    filt.add_argument("--total-sasa", type=float, default=80.0,
                      help="Min total SASA Å² (default: 80)")
    filt.add_argument("--apolar-frac", type=float, default=0.4,
                      help="Min apolar alpha-sphere fraction (default: 0.4)")
    filt.add_argument("--n-spheres", type=int, default=15,
                      help="Min alpha sphere count (default: 15)")

    ori = p.add_argument_group("orientation thresholds")
    ori.add_argument("--min-score", type=float, default=-5.0,
                     help="Max Vina score to accept kcal/mol (default: -5.0)")
    ori.add_argument("--max-ring-dist", type=float, default=8.0,
                     help="Max ring-to-Fe distance Å (default: 8.0)")
    ori.add_argument("--flexibility", type=float, default=2.0,
                     help="Orientation flexibility: d_ring < d_oxy + flexibility "
                          "(default: 2.0)")

    args = p.parse_args()

    input_dir   = Path(args.input_dir).resolve()
    outputs_dir = Path(args.outputs_dir).resolve()
    docking_out = Path(args.docking_out).resolve()
    ligand      = Path(args.ligand).resolve()

    if not input_dir.is_dir():
        sys.exit(f"input_dir not found: {input_dir}")
    if not ligand.exists():
        sys.exit(f"Ligand file not found: {ligand}")

    outputs_dir.mkdir(parents=True, exist_ok=True)
    docking_out.mkdir(parents=True, exist_ok=True)

    # ── Stage 0: receptor preparation ────────────────────────────────────
    if args.prepare:
        print("\n" + "=" * 60)
        print("STAGE 0 — RECEPTOR PREPARATION (PDB/CIF → PDBQT)")
        print("=" * 60)
        input_files = (sorted(input_dir.glob("*.pdb"))
                       + sorted(input_dir.glob("*.cif")))
        remove_water = not args.keep_water
        prepared = skipped_prep = 0
        for src in input_files:
            out_pdbqt = input_dir / src.with_suffix(".pdbqt").name
            if out_pdbqt.exists():
                print(f"  [SKIP] {src.name}  (PDBQT already exists)")
                skipped_prep += 1
            else:
                if convert(src, out_pdbqt, remove_water=remove_water):
                    prepared += 1
        print(f"\n{prepared} prepared, {skipped_prep} already existed.")
    else:
        missing = [
            f.stem for f in sorted(input_dir.glob("*.pdb")) + sorted(input_dir.glob("*.cif"))
            if not (input_dir / f.with_suffix(".pdbqt").name).exists()
        ]
        if missing:
            print(f"\n[INFO] {len(missing)} PDB(s) have no matching .pdbqt and will "
                  f"be skipped at the docking stage: {', '.join(missing)}")
            print("       Run with --prepare to auto-convert them.")

    # ── Stage 1: fpocket ──────────────────────────────────────────────────
    if not args.skip_fpocket:
        run(
            [sys.executable, SCRIPT_DIR / "run_fpocket.py",
             str(input_dir), "-o", str(outputs_dir)],
            "fpocket",
        )
    else:
        print("\n[fpocket] skipped — using existing outputs in", outputs_dir)

    # ── Stage 2: filter ───────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("STAGE 2 — POCKET FILTERING (heme-proximal)")
    print("=" * 60)

    rows = load_all_metrics(outputs_dir)  # used only to enumerate stems
    if not rows:
        sys.exit(f"No fpocket outputs found in {outputs_dir}")

    # all_results: (stem, status, pocket_n, metrics, fail_reason)
    all_results = []
    passing = []  # (stem, pocket_n, metrics)

    if args.skip_filter:
        for stem, m in rows:
            all_results.append((stem, "PASS", 1, m, ""))
            passing.append((stem, 1, m))
        print(f"Filter skipped — {len(passing)} structure(s) forwarded to docking "
              f"(defaulting to pocket 1).")
    else:
        for stem, _ in rows:
            try:
                receptor_pdb = find_receptor_pdb(stem, input_dir, outputs_dir)
                fe, _ = find_heme_fe(str(receptor_pdb))
            except (FileNotFoundError, ValueError) as e:
                reason = f"no heme: {e}"
                print(f"  FAIL  {stem:30s}  {reason}")
                all_results.append((stem, "FAIL", None, None, reason))
                continue

            pocket_n, metrics, fail_reason = find_heme_pocket(
                outputs_dir, stem, fe,
                fe_proximity=args.fe_proximity,
                druggability=args.druggability,
                volume=args.volume,
                total_sasa=args.total_sasa,
                apolar_frac=args.apolar_frac,
                n_spheres=args.n_spheres,
            )
            if pocket_n is None:
                print(f"  FAIL  {stem:30s}  {fail_reason}")
                all_results.append((stem, "FAIL", None, metrics, fail_reason))
            else:
                drugg = metrics.get("druggability_score", float("nan"))
                vol   = metrics.get("volume", float("nan"))
                print(f"  PASS  {stem:30s}  pocket {pocket_n}"
                      f"  drugg={drugg:.3f}  vol={vol:.0f}")
                all_results.append((stem, "PASS", pocket_n, metrics, ""))
                passing.append((stem, pocket_n, metrics))
        print(f"\n{len(passing)} / {len(rows)} structures pass.\n")

    # Build per-stem records that accumulate results from all stages
    def _fmt(m, key, fmt):
        if not m or key not in m:
            return ""
        try:
            return format(m[key], fmt)
        except (ValueError, TypeError):
            return ""

    records = {}
    for stem, status, pocket_n, metrics, fail_reason in all_results:
        records[stem] = {
            "filter_status": status,
            "pocket":        str(pocket_n) if pocket_n is not None else "",
            "druggability":  _fmt(metrics, "druggability_score", ".3f"),
            "volume":        _fmt(metrics, "volume",             ".1f"),
            "total_sasa":    _fmt(metrics, "total_sasa",         ".1f"),
            "apolar_frac":   _fmt(metrics, "apolar_sphere_frac", ".3f"),
            "n_spheres":     _fmt(metrics, "n_alpha_spheres",    ".0f"),
            "filter_reason": fail_reason,
            "docking_status": "not_reached" if status == "FAIL" else "",
            "n_poses":        "",
            "n_passing_poses":"",
            "best_score":     "",
            "best_d_ring":    "",
            "best_d_oxy":     "",
            "outcome":        "filtered_out" if status == "FAIL" else "",
        }

    def _write_tsv(path, records):
        cols = [
            "stem", "filter_status", "pocket",
            "druggability", "volume", "total_sasa", "apolar_frac", "n_spheres",
            "filter_reason",
            "docking_status", "n_poses", "n_passing_poses",
            "best_score", "best_d_ring", "best_d_oxy",
            "outcome",
        ]
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for stem, r in records.items():
                row = [stem] + [str(r.get(c, "")) for c in cols[1:]]
                fh.write("\t".join(row) + "\n")
        print(f"Results saved to {path}")

    if args.skip_docking:
        for stem in records:
            if records[stem]["docking_status"] == "":
                records[stem]["docking_status"] = "skipped"
                records[stem]["outcome"] = "skipped"
        _write_tsv(outputs_dir / "pocket_selection.tsv", records)
        print("\nDocking skipped (--skip-docking). Done.")
        sys.exit(0)

    # ── Stage 3: docking ──────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("STAGE 3 — DOCKING")
    print("=" * 60)

    docked = []   # (stem, out_pdbqt)
    n_skipped = 0

    for stem, pocket_n, _ in passing:
        receptor_pdbqt = input_dir / f"{stem}.pdbqt"
        if not receptor_pdbqt.exists():
            print(f"  [SKIP] {stem} — no matching .pdbqt in {input_dir}")
            records[stem]["docking_status"] = "no_pdbqt"
            records[stem]["outcome"] = "no_pdbqt"
            n_skipped += 1
            continue

        # Auto-derive box from the selected heme-proximal pocket
        try:
            center, size = derive_box(stem, input_dir, outputs_dir,
                                      args.buffer, args.heme_max_dist, pocket_n)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [SKIP] {stem} — box computation failed: {e}")
            records[stem]["docking_status"] = "box_failed"
            records[stem]["outcome"] = "box_failed"
            n_skipped += 1
            continue

        if args.box_size:
            size = args.box_size

        cx, cy, cz = [round(v, 2) for v in center]
        sx, sy, sz = [round(v, 1) for v in size]

        # Predict output filename before running (docking.py uses 1 decimal coords).
        coords_str = f"{cx:.1f}_{cy:.1f}_{cz:.1f}"
        rec_stem = receptor_pdbqt.stem.replace("_clean_h_meeko", "")
        lig_stem = ligand.stem
        out_pdbqt = docking_out / f"{rec_stem}_{lig_stem}_targeted_{coords_str}.pdbqt"
        if not out_pdbqt.exists():
            candidates = sorted(docking_out.glob(f"*{lig_stem}_targeted_{coords_str}*.pdbqt"))
            if candidates:
                out_pdbqt = candidates[-1]

        print(f"\n  {stem}")
        print(f"    pocket     : {pocket_n}")
        print(f"    box center : {cx} {cy} {cz}")
        print(f"    box size   : {sx} {sy} {sz}")

        if out_pdbqt.exists():
            print(f"    [SKIP] docking output already exists: {out_pdbqt.name}")
            records[stem]["docking_status"] = "already_docked"
            docked.append((stem, out_pdbqt))
            continue

        cmd = [
            sys.executable, SCRIPT_DIR / "docking.py",
            "-r", str(receptor_pdbqt),
            "-l", str(ligand),
            "-o", str(docking_out),
            "--center", cx, cy, cz,
            "--box_size", sx, sy, sz,
            "--exhaustiveness", args.exhaustiveness,
            "--n_poses", args.n_poses,
        ]
        try:
            run(cmd, f"vina/{stem}")
        except RuntimeError as e:
            print(f"  [ERROR] {e}")
            records[stem]["docking_status"] = "docking_failed"
            records[stem]["outcome"] = "docking_failed"
            n_skipped += 1
            continue

        if out_pdbqt.exists():
            records[stem]["docking_status"] = "docked"
            docked.append((stem, out_pdbqt))
        else:
            print(f"  [WARN] docking output not found for {stem}")
            records[stem]["docking_status"] = "docking_failed"
            records[stem]["outcome"] = "docking_failed"
            n_skipped += 1

    # ── Stage 4: pose analysis ────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("STAGE 4 — POSE ORIENTATION ANALYSIS")
    print("=" * 60)

    summary_rows = []

    for stem, pdbqt_path in docked:
        try:
            receptor_pdb = find_receptor_pdb(stem, input_dir, outputs_dir)
            fe, ref = find_heme_fe(str(receptor_pdb))
        except (FileNotFoundError, ValueError) as e:
            print(f"  [SKIP] {stem} — cannot locate heme for orientation check: {e}")
            records[stem]["outcome"] = "heme_not_found"
            continue

        results = analyze_file(
            str(pdbqt_path), fe,
            args.min_score, args.max_ring_dist, args.flexibility,
            top_n=None,
        )
        print_results(pdbqt_path, results, args.min_score,
                      args.max_ring_dist, args.flexibility)

        n_pass = sum(r["pass"] for r in results)
        best_pass = next((r for r in results if r["pass"]), None)

        records[stem]["n_poses"]         = str(len(results))
        records[stem]["n_passing_poses"] = str(n_pass)
        if best_pass:
            records[stem]["best_score"]  = f"{best_pass['score']:.2f}"
            records[stem]["best_d_ring"] = f"{best_pass['d_ring']:.2f}"
            d_oxy = best_pass.get("d_oxy")
            records[stem]["best_d_oxy"]  = f"{d_oxy:.2f}" if d_oxy is not None else ""
            records[stem]["outcome"]     = "pass"
        else:
            records[stem]["outcome"]     = "no_passing_poses"

        summary_rows.append({
            "stem": stem,
            "n_poses": len(results),
            "n_pass": n_pass,
            "best_pass": best_pass,
        })

    # ── Write final TSV ───────────────────────────────────────────────────
    _write_tsv(outputs_dir / "pocket_selection.tsv", records)

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    print(f"  Input structures   : {len(rows)}")
    print(f"  After filtering    : {len(passing)}")
    print(f"  Docking completed  : {len(docked)}")
    print(f"  Docking skipped    : {n_skipped}")
    n_with_pass = sum(1 for r in summary_rows if r["n_pass"] > 0)
    print(f"  Orientation pass   : {n_with_pass} structure(s) with ≥1 passing pose\n")

    if summary_rows:
        print(f"  {'Structure':<30} {'Poses':>6} {'Pass':>5}  Best passing pose")
        print(f"  {'─'*70}")
        for row in summary_rows:
            bp = row["best_pass"]
            if bp:
                d_oxy = bp.get("d_oxy")
                detail = (f"pose {bp['pose']}  score={bp['score']:.2f}  "
                          f"d_ring={bp['d_ring']:.1f}Å"
                          + (f"  d_oxy={d_oxy:.1f}Å" if d_oxy is not None else ""))
            else:
                detail = "—"
            print(f"  {row['stem']:<30} {row['n_poses']:>6} {row['n_pass']:>5}  {detail}")


if __name__ == "__main__":
    main()
