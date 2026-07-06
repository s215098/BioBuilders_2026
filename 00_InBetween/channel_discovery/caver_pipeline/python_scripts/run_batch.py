import argparse
import csv
import math
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from extract_tunnel_R import get_tunnel_axis_values, plot_distribution, plot_tunnel_shape

ROOT = Path(__file__).parent.parent
CAVER_JAR  = ROOT / 'caver' / 'caver.jar'
CAVER_HOME = ROOT / 'caver'
CAVER_LIB  = ROOT / 'caver' / 'lib'

HEME_NAMES = {'HEM', 'HEME', 'HEC', 'HEA', 'HEB', 'HEG', 'HDD', 'HAS', 'LIG'}


# ── HEME / starting-point detection ──────────────────────────────────────────

def find_heme(pdb_file):
    """
    Return (fe_coords, protein_chain) for the heme iron, or (None, None) if not found.
    fe_coords is (x, y, z) of the Fe atom; protein_chain is the chain containing ATOM records.
    """
    fe_coords = None
    protein_chain = None
    max_res = 0

    with open(pdb_file) as f:
        for line in f:
            if line.startswith('HETATM'):
                res_name = line[17:21].strip()
                if res_name in HEME_NAMES and line[12:16].strip() in ('FE', 'FE2', 'FE3'):
                    fe_coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            elif line.startswith('ATOM'):
                protein_chain = line[21]
                try:
                    max_res = max(max_res, int(line[22:26]))
                except ValueError:
                    pass

    if fe_coords is None or protein_chain is None:
        return None, None

    return fe_coords, protein_chain, max_res


def write_pdb_with_fe_atom(pdb_file, out_path, fe_coords, protein_chain, max_res):
    """
    Copy pdb_file to out_path, appending a synthetic ATOM record for the Fe so CAVER
    can resolve it as starting_point_residue. Returns (res_num_str, chain).
    """
    fe_res = max_res + 1
    fx, fy, fz = fe_coords
    fe_line = (f"ATOM  {9999:5d}  FE  FEX {protein_chain}{fe_res:4d}    "
               f"{fx:8.3f}{fy:8.3f}{fz:8.3f}  1.00  0.00          FE  \n")

    if out_path.is_symlink() or out_path.exists():
        out_path.unlink()

    lines = open(pdb_file).readlines()
    with open(out_path, 'w') as f:
        for line in lines:
            stripped = line.strip()
            if stripped == 'END' or stripped == '':
                continue
            f.write(line)
        f.write(fe_line)
        f.write('END\n')

    return str(fe_res), protein_chain


# ── fpocket integration ───────────────────────────────────────────────────────

def load_pocket_selection(fpocket_dir):
    """
    Parse pocket_selection.tsv and return {stem: pocket_number} for PASS entries only.
    """
    tsv = Path(fpocket_dir) / 'pocket_selection.tsv'
    if not tsv.exists():
        return {}

    selection = {}
    with open(tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            stem   = row.get('stem', '').strip()
            status = row.get('filter_status', '').strip()
            pocket = row.get('pocket', '').strip()
            if status == 'PASS' and stem and pocket:
                try:
                    selection[stem] = int(pocket)
                except ValueError:
                    pass
    return selection


def pocket_center(fpocket_dir, stem, pocket_num):
    """
    Return (x, y, z) center of mass of alpha spheres for the given pocket,
    or None if the file is not found.
    """
    pqr = Path(fpocket_dir) / f'{stem}_out' / 'pockets' / f'pocket{pocket_num}_vert.pqr'
    if not pqr.exists():
        return None

    coords = []
    with open(pqr) as f:
        for line in f:
            if line.startswith('ATOM'):
                fields = line.split()
                coords.append((float(fields[5]), float(fields[6]), float(fields[7])))

    if not coords:
        return None

    n = len(coords)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)


def tunnel_pocket_dist(csv_path, tunnel_id, pocket_ctr):
    """
    Return the minimum distance (Å) between any point on the tunnel's 3D path
    and pocket_ctr.  Returns None if XYZ data is missing.
    """
    x_vals = get_tunnel_axis_values(str(csv_path), tunnel_id, 'X')
    y_vals = get_tunnel_axis_values(str(csv_path), tunnel_id, 'Y')
    z_vals = get_tunnel_axis_values(str(csv_path), tunnel_id, 'Z')

    if not x_vals or not y_vals or not z_vals:
        return None

    px, py, pz = pocket_ctr
    n = min(len(x_vals), len(y_vals), len(z_vals))
    dists = [
        math.sqrt((x_vals[i] - px) ** 2 +
                  (y_vals[i] - py) ** 2 +
                  (z_vals[i] - pz) ** 2)
        for i in range(n)
    ]
    return sum(dists) / len(dists)


def read_tunnel_characteristics(caver_out_dir):
    """
    Parse tunnel_characteristics.csv and return
    {tunnel_id: {bottleneck_radius, length, throughput, curvature}}.
    """
    path = caver_out_dir / 'analysis' / 'tunnel_characteristics.csv'
    if not path.exists():
        return {}

    chars = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for row in reader:
            try:
                tid = int(row['Tunnel'].strip())
            except (KeyError, ValueError):
                continue
            chars[tid] = {
                'bottleneck_radius': row.get('Bottleneck radius', '').strip(),
                'length':            row.get('Length', '').strip(),
                'throughput':        row.get('Throughput', '').strip(),
                'curvature':         row.get('Curvature', '').strip(),
            }
    return chars


# ── config generation ─────────────────────────────────────────────────────────

def make_config(template_path, out_path, res_num, chain, probe_radius=None):
    """Write a per-PDB config, overriding starting_point_residue and optionally probe_radius."""
    with open(template_path) as f:
        lines = f.readlines()

    replaced_spr = False
    replaced_pr  = False
    out_lines = []
    for line in lines:
        if re.match(r'\s*starting_point_residue\b', line):
            out_lines.append(f'starting_point_residue     {res_num} {chain}\n')
            replaced_spr = True
        elif probe_radius is not None and re.match(r'\s*probe_radius\b', line):
            out_lines.append(f'probe_radius                {probe_radius}\n')
            replaced_pr = True
        else:
            out_lines.append(line)

    if not replaced_spr:
        out_lines.append(f'\nstarting_point_residue     {res_num} {chain}\n')
    if probe_radius is not None and not replaced_pr:
        out_lines.append(f'\nprobe_radius                {probe_radius}\n')

    with open(out_path, 'w') as f:
        f.writelines(out_lines)


# ── CAVER execution ───────────────────────────────────────────────────────────

def run_caver(pdb_file, pdb_input_dir, config, out_dir):
    """Run CAVER. pdb_input_dir must contain exactly one PDB file (real file or symlink)."""
    pdb_input_dir.mkdir(parents=True, exist_ok=True)
    target = pdb_input_dir / pdb_file.name
    # Only symlink if not already a real file (modified PDB written there directly)
    if not (target.exists() and not target.is_symlink()):
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(pdb_file.resolve())

    cmd = [
        'java', '-Xmx1200m',
        '-cp', str(CAVER_LIB),
        '-jar', str(CAVER_JAR),
        '-home', str(CAVER_HOME),
        '-pdb', str(pdb_input_dir),
        '-conf', str(config),
        '-out', str(out_dir),
    ]
    print('  Running CAVER...', flush=True)
    subprocess.run(cmd, capture_output=True, text=True)

    log = out_dir / 'log.txt'
    if log.exists() and ('Finished successfully' in log.read_text() or
                         'Calculation finished' in log.read_text()):
        return True

    print('  ERROR: CAVER did not finish successfully — check log:')
    print(f'    {log}')
    return False


# ── analysis & plotting ───────────────────────────────────────────────────────

def tunnel_ids_from_csv(csv_path):
    ids = set()
    with open(csv_path, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            row = [c.strip() for c in row]
            if len(row) < 13:
                continue
            try:
                ids.add(int(row[2]))
            except ValueError:
                pass
    return sorted(ids)


def run_analysis(caver_out_dir, plots_dir, ylim, pocket_ctr=None, pocket_radius=12.0,
                 max_curvature=None, max_length=None, make_plots=True):
    """
    Evaluate the pocket-distance / curvature / length filters for every tunnel and,
    when make_plots is True, save plots for the kept ones. Returns a list of row dicts
    for tunnel_selection.csv, sorted kept-first then by throughput descending. With
    make_plots=False the rows are still computed but no plot files are written (used to
    add an already-computed structure to the CSV without regenerating its plots).
    """
    csv_path = caver_out_dir / 'analysis' / 'tunnel_profiles.csv'
    if not csv_path.exists():
        print('  WARNING: tunnel_profiles.csv not found — CAVER may have found no tunnels')
        return []

    if make_plots:
        plots_dir.mkdir(parents=True, exist_ok=True)
    chars   = read_tunnel_characteristics(caver_out_dir)
    all_ids = tunnel_ids_from_csv(csv_path)
    kept, skipped = [], []
    rows = []

    for tid in all_ids:
        c = chars.get(tid, {})
        reasons = []

        # ── pocket distance filter ────────────────────────────────────────────
        mean_dist = None
        if pocket_ctr is not None:
            mean_dist = tunnel_pocket_dist(csv_path, tid, pocket_ctr)
            if mean_dist is not None and mean_dist > pocket_radius:
                reasons.append(f'mean_dist {mean_dist:.1f} Å > {pocket_radius} Å')

        # ── curvature filter ──────────────────────────────────────────────────
        curv = None
        if max_curvature is not None and c.get('curvature'):
            try:
                curv = float(c['curvature'])
                if curv > max_curvature:
                    reasons.append(f'curvature {curv:.2f} > {max_curvature}')
            except ValueError:
                pass

        # ── length filter ─────────────────────────────────────────────────────
        if max_length is not None and c.get('length'):
            try:
                length = float(c['length'])
                if length > max_length:
                    reasons.append(f'length {length:.1f} Å > {max_length} Å')
            except ValueError:
                pass

        status = 'discarded' if reasons else 'kept'
        (skipped if reasons else kept).append(tid)
        rows.append({
            'tunnel':        tid,
            'filter_status': status,
            'discard_reason': '; '.join(reasons),
            'pocket_dist':   f'{mean_dist:.2f}' if mean_dist is not None else '',
            **c,
        })

        if not reasons and make_plots:
            r_values    = get_tunnel_axis_values(str(csv_path), tid, 'R')
            dist_values = get_tunnel_axis_values(str(csv_path), tid, 'distance')
            if r_values:
                plot_distribution(r_values, tid, str(plots_dir))
            if r_values and dist_values:
                plot_tunnel_shape(r_values, dist_values, tid, str(plots_dir), ylim=ylim)

    if pocket_ctr is not None or max_curvature is not None:
        kept_str    = ', '.join(str(t) for t in kept)    or 'none'
        skipped_str = ', '.join(str(t) for t in skipped) or 'none'
        print(f'  Filter: kept [{kept_str}], discarded [{skipped_str}]')

    print(f'  Plots saved to {plots_dir}' if make_plots else '  Plots skipped (already computed)')

    # Sort: kept first, then by throughput descending within each group
    def sort_key(row):
        order = 0 if row['filter_status'] == 'kept' else 1
        try:
            tp = -float(row.get('throughput', 0))
        except (ValueError, TypeError):
            tp = 0
        return (order, tp)

    return sorted(rows, key=sort_key)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Batch-run CAVER on every PDB in a folder, then generate plots for each.'
    )
    parser.add_argument('--pdb',           type=Path,  default=None,
                        help='Run on a single PDB file instead of a whole folder')
    parser.add_argument('--pdb-dir',       type=Path,  default=ROOT / 'input_files',
                        help='Folder containing .pdb files (default: input_files/)')
    parser.add_argument('--config',        type=Path,  default=ROOT / 'config.txt',
                        help='Config template (default: config.txt)')
    parser.add_argument('--results-dir',   type=Path,  default=ROOT / 'results',
                        help='Root folder for all outputs (default: results/)')
    parser.add_argument('--fpocket-dir',   type=Path,  default=None,
                        help='fpocket outputs folder containing pocket_selection.tsv')
    parser.add_argument('--pocket-radius',  type=float, default=12.0,
                        help='Max mean distance (Å) of tunnel path from pocket center (default: 12)')
    parser.add_argument('--max-curvature', type=float, default=1.5,
                        help='Discard tunnels with curvature above this value (default: 1.5)')
    parser.add_argument('--max-length',   type=float, default=25.0,
                        help='Discard tunnels longer than this value in Å (default: 25)')
    parser.add_argument('--probe-radius',  type=float, default=None,
                        help='Override probe radius in config (default: use config.txt value)')
    parser.add_argument('--ylim',          type=float, default=5,
                        help='Y-axis half-range for tunnel shape plots (default: 5)')
    parser.add_argument('--no-plots',      action='store_true',
                        help='Skip plot generation after CAVER runs')
    parser.add_argument('--no-caver',      action='store_true',
                        help='Skip CAVER, only (re-)generate plots from existing output')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip both the CAVER run and plot regeneration for structures that '
                             'already have output (results/<stem>/out/summary.txt); their rows are '
                             'still added to the combined CSV, so new PDBs are computed and added '
                             'without recomputing or re-plotting the existing ones')
    args = parser.parse_args()

    # Load fpocket pocket selection if provided
    pocket_selection = {}
    if args.fpocket_dir:
        pocket_selection = load_pocket_selection(args.fpocket_dir)
        print(f'fpocket: {len(pocket_selection)} PASS structure(s): '
              f'{", ".join(f"{s}→pocket{n}" for s, n in pocket_selection.items())}')

    # Collect PDB files
    if args.pdb:
        if not args.pdb.is_file():
            print(f'File not found: {args.pdb}')
            sys.exit(1)
        pdbs = [args.pdb]
    else:
        pdbs = sorted(args.pdb_dir.glob('*.pdb'))
        if not pdbs:
            print(f'No .pdb files found in {args.pdb_dir}')
            sys.exit(1)

    print(f'Found {len(pdbs)} PDB file(s)')

    all_tunnel_rows = []

    for pdb in pdbs:
        stem          = pdb.stem
        pdb_dir       = args.results_dir / stem
        pdb_input_dir = pdb_dir / 'pdb'
        caver_out     = pdb_dir / 'out'
        plots_out     = pdb_dir / 'plots'
        config        = pdb_dir / 'config.txt'

        print(f'\n[{stem}]')

        if args.fpocket_dir and stem not in pocket_selection:
            print(f'  Skipped — not in fpocket PASS list')
            continue

        caver_out.mkdir(parents=True, exist_ok=True)

        caver_skipped = False
        if not args.no_caver:
            if args.skip_existing and (caver_out / 'summary.txt').exists():
                print('  Skipping CAVER — output already exists (--skip-existing)')
                caver_skipped = True
            else:
                result = find_heme(pdb)
                if result[0] is None:
                    print(f'  WARNING: no HEME Fe found in {pdb.name}, falling back to template config')
                    config = args.config
                    ok = run_caver(pdb, pdb_input_dir, config, caver_out)
                else:
                    fe_coords, protein_chain, max_res = result
                    print(f'  Fe at ({fe_coords[0]:.2f}, {fe_coords[1]:.2f}, {fe_coords[2]:.2f})')
                    pdb_input_dir.mkdir(parents=True, exist_ok=True)
                    modified_pdb = pdb_input_dir / pdb.name
                    fe_res_num, fe_chain = write_pdb_with_fe_atom(
                        pdb, modified_pdb, fe_coords, protein_chain, max_res)
                    print(f'  Starting point: synthetic Fe residue {fe_res_num} {fe_chain}')
                    make_config(args.config, config, fe_res_num, fe_chain,
                                probe_radius=args.probe_radius)
                    ok = run_caver(modified_pdb, pdb_input_dir, config, caver_out)

                if not ok:
                    print('  Skipping plots due to CAVER error')
                    continue

        if not args.no_plots:
            # Resolve pocket center for this structure if available
            pctr = None
            pnum = None
            if args.fpocket_dir and stem in pocket_selection:
                pnum = pocket_selection[stem]
                pctr = pocket_center(args.fpocket_dir, stem, pnum)
                if pctr:
                    print(f'  Pocket filter: using fpocket pocket {pnum} '
                          f'(center {pctr[0]:.1f}, {pctr[1]:.1f}, {pctr[2]:.1f}), '
                          f'radius {args.pocket_radius} Å')
                else:
                    print(f'  WARNING: pocket {pnum} PQR not found for {stem}, no filtering')
            elif args.fpocket_dir:
                print(f'  No fpocket PASS entry for {stem} — keeping all tunnels')

            rows = run_analysis(caver_out, plots_out, args.ylim,
                                pocket_ctr=pctr, pocket_radius=args.pocket_radius,
                                max_curvature=args.max_curvature,
                                max_length=args.max_length,
                                make_plots=not caver_skipped)
            fpocket_status = 'PASS' if (args.fpocket_dir and stem in pocket_selection) else ''
            for row in rows:
                all_tunnel_rows.append({
                    'stem':           stem,
                    'fpocket_status': fpocket_status,
                    'pocket_num':     pnum if pnum is not None else '',
                    **row,
                })

    print('\nDone.')

    if all_tunnel_rows:
        csv_out = args.results_dir / 'tunnel_selection.csv'
        fieldnames = ['stem', 'fpocket_status', 'pocket_num', 'tunnel',
                      'filter_status', 'discard_reason', 'throughput',
                      'bottleneck_radius', 'length', 'curvature', 'pocket_dist']
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_out, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_tunnel_rows)
        print(f'Tunnel selection summary → {csv_out}')


if __name__ == '__main__':
    main()
