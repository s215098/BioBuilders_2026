import csv
import argparse
import os
import matplotlib.pyplot as plt

CSV_DEFAULT = os.path.join(os.path.dirname(__file__), '..', 'out', 'analysis', 'tunnel_profiles.csv')


def get_tunnel_axis_values(csv_file, tunnel_id, axis):
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header

        for row in reader:
            row = [cell.strip() for cell in row]

            if len(row) < 13:
                continue

            try:
                row_tunnel = int(row[2])
            except ValueError:
                continue

            if row_tunnel == tunnel_id and row[12] == axis:
                return [float(v) for v in row[13:] if v and v != '-']

    return None


def plot_distribution(values, tunnel_id, out_dir):
    fig, ax = plt.subplots()
    ax.hist(values, bins='auto', edgecolor='black')
    ax.set_xlabel('R (Å)')
    ax.set_ylabel('Count')
    ax.set_title(f'Tunnel {tunnel_id} — R value distribution')
    path = os.path.join(out_dir, f'tunnel_{tunnel_id}_R_distribution.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {path}')


def plot_tunnel_shape(r_values, dist_values, tunnel_id, out_dir, ylim=3):
    # Round distances to 2 decimals for x-axis labels
    dist_rounded = [round(d, 2) for d in dist_values]

    n = min(len(r_values), len(dist_values))
    x = dist_rounded[:n]
    r = r_values[:n]

    upper = r
    lower = [-v for v in r]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, upper, color='steelblue', linewidth=1.5)
    ax.plot(x, lower, color='steelblue', linewidth=1.5)
    ax.fill_between(x, lower, upper, alpha=0.25, color='steelblue')

    ax.set_ylim(-ylim, ylim)
    ax.set_xlabel('Distance along tunnel (Å)')
    ax.set_ylabel('R (Å)')
    ax.set_title(f'Tunnel {tunnel_id} — cross-section profile')
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')

    path = os.path.join(out_dir, f'tunnel_{tunnel_id}_R_shape.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract R values for a given tunnel from tunnel_profiles.csv')
    parser.add_argument('--tunnel', type=int, default=1, help='Tunnel ID (default: 1)')
    parser.add_argument('--csv', type=str, default=CSV_DEFAULT, help='Path to tunnel_profiles.csv')
    parser.add_argument('--out', type=str, default=os.path.dirname(__file__), help='Output directory for plots')
    parser.add_argument('--ylim', type=float, default=3, help='Y-axis range: plot goes from -ylim to +ylim (default: 3)')
    args = parser.parse_args()

    r_values = get_tunnel_axis_values(args.csv, args.tunnel, 'R')
    dist_values = get_tunnel_axis_values(args.csv, args.tunnel, 'distance')

    if r_values is None:
        print(f"Tunnel {args.tunnel} with Axis 'R' not found in {args.csv}")
    else:
        print(f"R values for tunnel {args.tunnel} ({len(r_values)} points):")
        print(r_values)

        plot_distribution(r_values, args.tunnel, args.out)

        if dist_values is None:
            print(f"Warning: distance row not found for tunnel {args.tunnel}, skipping shape plot")
        else:
            print(f"\nDistance values for tunnel {args.tunnel} ({len(dist_values)} points):")
            print([round(d, 2) for d in dist_values])
            plot_tunnel_shape(r_values, dist_values, args.tunnel, args.out, ylim=args.ylim)
