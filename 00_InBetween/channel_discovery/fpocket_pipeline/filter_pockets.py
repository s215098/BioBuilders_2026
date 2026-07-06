#!/usr/bin/env python3
"""Filter structures by fpocket Pocket 1 metrics.

Scans all <stem>_out/ directories under outputs_dir, reads the per-pocket
scores from <stem>_info.txt, and reports which structures pass the threshold
filters on Pocket 1.

Defaults are intentionally loose — the purpose is to discard obvious non-hits
(tiny or buried sites) before docking, not to be selective.

Examples
--------
    python3 filter_pockets.py
    python3 filter_pockets.py outputs/ --druggability 0.3 --volume 400
    python3 filter_pockets.py --out-tsv metrics.tsv --out-list passing.txt
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np

# Allow importing read_coords from pocket_box when used as a library
_HERE = Path(__file__).parent.resolve()
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from pocket_box import read_coords

FIELD_MAP = {
    "Score":                        "score",
    "Druggability Score":           "druggability_score",
    "Number of Alpha Spheres":      "n_alpha_spheres",
    "Total SASA":                   "total_sasa",
    "Polar SASA":                   "polar_sasa",
    "Apolar SASA":                  "apolar_sasa",
    "Volume":                       "volume",
    "Apolar alpha sphere proportion": "apolar_sphere_frac",
    "Hydrophobicity score":         "hydrophobicity_score",
    "Polarity score":               "polarity_score",
    "Alpha sphere density":         "alpha_sphere_density",
    "Flexibility":                  "flexibility",
}

DISPLAY_COLS = [
    ("score",             "Score"),
    ("druggability_score","Drugg"),
    ("volume",            "Vol(Å³)"),
    ("total_sasa",        "SASA(Å²)"),
    ("apolar_sasa",       "ApolSASA"),
    ("apolar_sphere_frac","ApolFrac"),
    ("n_alpha_spheres",   "nSph"),
]


def parse_pocket_info(path):
    """Parse fpocket *_info.txt. Returns {pocket_num: {field: value}}."""
    pockets = {}
    current_n = None
    current = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^Pocket (\d+) :", line)
            if m:
                if current_n is not None:
                    pockets[current_n] = current
                current_n = int(m.group(1))
                current = {}
                continue
            if ":" in line and current_n is not None:
                key, _, val = line.partition(":")
                key = key.strip()
                mapped = FIELD_MAP.get(key)
                if mapped is None:
                    continue
                try:
                    current[mapped] = float(val.strip())
                except ValueError:
                    pass
    if current_n is not None:
        pockets[current_n] = current
    return pockets


def load_all_metrics(outputs_dir):
    """Return list of (stem, pocket1_metrics) for every _out/ folder found."""
    results = []
    for out_dir in sorted(Path(outputs_dir).iterdir()):
        if not out_dir.is_dir() or not out_dir.name.endswith("_out"):
            continue
        stem = out_dir.name[:-4]
        info_file = out_dir / f"{stem}_info.txt"
        if not info_file.exists():
            print(f"[WARN] no info file: {info_file}", file=sys.stderr)
            continue
        pockets = parse_pocket_info(info_file)
        if 1 not in pockets:
            print(f"[WARN] no Pocket 1 in {info_file}", file=sys.stderr)
            continue
        results.append((stem, pockets[1]))
    return results


def apply_filters(metrics, druggability, volume, total_sasa,
                  apolar_frac, n_spheres):
    reasons = []
    if metrics.get("druggability_score", 0) < druggability:
        reasons.append(f"druggability {metrics.get('druggability_score', '?'):.3f} < {druggability}")
    if metrics.get("volume", 0) < volume:
        reasons.append(f"volume {metrics.get('volume', '?'):.1f} < {volume}")
    if metrics.get("total_sasa", 0) < total_sasa:
        reasons.append(f"total_sasa {metrics.get('total_sasa', '?'):.1f} < {total_sasa}")
    if metrics.get("apolar_sphere_frac", 0) < apolar_frac:
        reasons.append(f"apolar_frac {metrics.get('apolar_sphere_frac', '?'):.2f} < {apolar_frac}")
    if metrics.get("n_alpha_spheres", 0) < n_spheres:
        reasons.append(f"n_spheres {int(metrics.get('n_alpha_spheres', 0))} < {n_spheres}")
    return reasons  # empty = passes


def find_heme_pocket(outputs_dir, stem, fe_coord,
                     fe_proximity=10.0,
                     druggability=0.2, volume=300.0, total_sasa=80.0,
                     apolar_frac=0.4, n_spheres=15,
                     max_pockets=20):
    """Return (pocket_n, metrics, fail_reason) for the first pocket that is
    heme-proximal and passes all metric filters.

    On success: (pocket_n, metrics, None)
    On failure: (None, best_candidate_metrics_or_None, reason_string)

    Iterates pockets in fpocket order (1, 2, 3, ...) and accepts a pocket when:
      1. Min distance from any alpha sphere to fe_coord <= fe_proximity
      2. total_sasa > 0  (has some solvent exposure — not fully buried)
      3. Passes volume / druggability / apolar_frac / n_spheres thresholds
    """
    out_dir = Path(outputs_dir) / f"{stem}_out"
    pockets_dir = out_dir / "pockets"
    info_file = out_dir / f"{stem}_info.txt"

    if not info_file.exists() or not pockets_dir.is_dir():
        return None, None, "no fpocket output"

    pockets_info = parse_pocket_info(info_file)

    pqr_files = sorted(pockets_dir.glob("pocket*_vert.pqr"))
    pocket_nums = []
    for pf in pqr_files:
        m = re.match(r"pocket(\d+)_vert\.pqr", pf.name)
        if m:
            pocket_nums.append(int(m.group(1)))
    pocket_nums.sort()

    fe = np.asarray(fe_coord)
    best_proximal = None  # (pocket_n, min_dist, metrics) — closest to Fe even if failing metrics

    for pocket_n in pocket_nums[:max_pockets]:
        pqr = pockets_dir / f"pocket{pocket_n}_vert.pqr"
        spheres = read_coords(str(pqr))
        if not len(spheres):
            continue

        min_dist = float(np.linalg.norm(spheres - fe, axis=1).min())
        if min_dist > fe_proximity:
            continue

        metrics = pockets_info.get(pocket_n, {})
        if metrics.get("total_sasa", 0) <= 0:
            continue

        if best_proximal is None:
            best_proximal = (pocket_n, min_dist, metrics)

        reasons = apply_filters(metrics, druggability, volume, total_sasa,
                                apolar_frac, n_spheres)
        if not reasons:
            return pocket_n, metrics, None

    # No pocket passed — build a descriptive reason
    if best_proximal is None:
        return None, None, f"no pocket within {fe_proximity}Å of Fe"
    pocket_n, min_dist, metrics = best_proximal
    reasons = apply_filters(metrics, druggability, volume, total_sasa,
                            apolar_frac, n_spheres)
    return None, metrics, (f"pocket {pocket_n} (min_dist={min_dist:.1f}Å): "
                           f"{'; '.join(reasons)}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("outputs_dir", nargs="?", default="outputs",
                   help="Directory containing <stem>_out/ folders (default: outputs/)")

    g = p.add_argument_group("filter thresholds (Pocket 1)")
    g.add_argument("--druggability", type=float, default=0.2,
                   help="Min druggability score (default: 0.2)")
    g.add_argument("--volume", type=float, default=500.0,
                   help="Min pocket volume in Å³ (default: 500)")
    g.add_argument("--total-sasa", type=float, default=80.0,
                   help="Min total SASA in Å² (default: 80)")
    g.add_argument("--apolar-frac", type=float, default=0.4,
                   help="Min apolar alpha-sphere fraction (default: 0.4)")
    g.add_argument("--n-spheres", type=int, default=15,
                   help="Min number of alpha spheres (default: 15)")

    p.add_argument("--out-tsv", metavar="FILE",
                   help="Write full metrics table to TSV")
    p.add_argument("--out-list", metavar="FILE",
                   help="Write passing stem names (one per line) to file")
    args = p.parse_args()

    outputs_dir = Path(args.outputs_dir)
    if not outputs_dir.is_dir():
        sys.exit(f"Outputs directory not found: {outputs_dir}")

    rows = load_all_metrics(outputs_dir)
    if not rows:
        sys.exit(f"No _out/ directories with info files found in {outputs_dir}")

    passing = []
    failing = []
    for stem, m in rows:
        reasons = apply_filters(
            m, args.druggability, args.volume, args.total_sasa,
            args.apolar_frac, args.n_spheres,
        )
        if reasons:
            failing.append((stem, m, reasons))
        else:
            passing.append((stem, m))

    # --- print table -----------------------------------------------------------
    header = f"{'Structure':<30}" + "".join(f"{lbl:>10}" for _, lbl in DISPLAY_COLS) + "  Status"
    sep = "-" * len(header)
    print(sep)
    print(header)
    print(sep)
    for stem, m in passing:
        vals = "".join(f"{m.get(k, float('nan')):10.3f}" for k, _ in DISPLAY_COLS)
        print(f"{stem:<30}{vals}  PASS")
    for stem, m, reasons in failing:
        vals = "".join(f"{m.get(k, float('nan')):10.3f}" for k, _ in DISPLAY_COLS)
        print(f"{stem:<30}{vals}  FAIL  ({'; '.join(reasons)})")
    print(sep)
    print(f"{len(passing)} / {len(rows)} structures pass.\n")

    # --- optional TSV ----------------------------------------------------------
    if args.out_tsv:
        tsv_cols = [k for k, _ in DISPLAY_COLS]
        with open(args.out_tsv, "w") as fh:
            fh.write("stem\t" + "\t".join(tsv_cols) + "\tpass\n")
            for stem, m in passing:
                vals = "\t".join(f"{m.get(k, '')}" for k in tsv_cols)
                fh.write(f"{stem}\t{vals}\t1\n")
            for stem, m, _ in failing:
                vals = "\t".join(f"{m.get(k, '')}" for k in tsv_cols)
                fh.write(f"{stem}\t{vals}\t0\n")
        print(f"Metrics table written to {args.out_tsv}")

    # --- optional passing list -------------------------------------------------
    if args.out_list:
        with open(args.out_list, "w") as fh:
            for stem, _ in passing:
                fh.write(stem + "\n")
        print(f"Passing stems written to {args.out_list}")

    return passing  # usable when imported


if __name__ == "__main__":
    main()
