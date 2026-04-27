import argparse
import os
import subprocess
import tempfile
import numpy as np
from pathlib import Path
from vina import Vina
from Bio.PDB import MMCIFParser, PDBIO


def prepare_receptor_mgltools(input_path, mgltools_python, prepare_receptor_script):
    """
    Converts CIF/PDB → PDBQT using MGLTools.
    Returns path to temporary PDBQT file.
    """
    input_path = Path(input_path)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)

        # Step 1: Convert CIF → PDB if needed
        if input_path.suffix.lower() == ".cif":
            pdb_path = tmp_dir / (input_path.stem + ".pdb")

            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("rec", str(input_path))

            io = PDBIO()
            io.set_structure(structure)
            io.save(str(pdb_path))

        else:
            pdb_path = input_path

        # Step 2: Convert PDB → PDBQT
        out_pdbqt = tmp_dir / (pdb_path.stem + ".pdbqt")

        cmd = [
            mgltools_python,
            prepare_receptor_script,
            "-r", str(pdb_path),
            "-o", str(out_pdbqt),
            "-A", "hydrogens",
        ]

        print(f"Preparing receptor: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"prepare_receptor4.py failed:\n{result.stderr}"
            )

        if not out_pdbqt.exists():
            raise FileNotFoundError("Receptor PDBQT not created.")

        # Copy to persistent temp file (since tmp_dir will vanish)
        final_path = tmp_dir / ("final_receptor.pdbqt")
        final_path.write_text(out_pdbqt.read_text())

        return final_path
    


def get_receptor_center_and_size(pdb_file):
    """Calculates the center and bounding box size of the receptor for blind docking."""
    coordinates = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinates.append([x, y, z])

    coords = np.array(coordinates)
    center = np.mean(coords, axis=0)

    # (Max - Min) per axis gives the full span of the protein
    size = np.max(coords, axis=0) - np.min(coords, axis=0)

    # Add a 10 Angstrom buffer
    size = size + 10.0

    return center.tolist(), size.tolist()


def prepare_ligand_mgltools(sdf_path, mgltools_python, prepare_ligand_script):
    """
    Converts a ligand SDF file to PDBQT format using MGLTools' prepare_ligand4.py.

    MGLTools runs under Python 2 and must be called as a subprocess from this
    Python 3 script. The resulting PDBQT content is written to a temporary file
    and read back as a string for use with Vina.

    Args:
        sdf_path (str | Path): Path to the input ligand SDF file.
        mgltools_python (str): Path to the MGLTools Python 2 interpreter,
            e.g. /opt/mgltools/bin/pythonsh
        prepare_ligand_script (str): Path to prepare_ligand4.py,
            e.g. /opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py

    Returns:
        str: The PDBQT file contents as a string, ready to pass to
             Vina.set_ligand_from_string().

    Raises:
        RuntimeError: If MGLTools exits with a non-zero return code.
        FileNotFoundError: If the expected PDBQT output file is not produced.
    """
    sdf_path = Path(sdf_path)

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_pdbqt = Path(tmp_dir) / (sdf_path.stem + ".pdbqt")

        cmd = [
            mgltools_python,
            prepare_ligand_script,
            "-l", str(sdf_path),
            "-o", str(out_pdbqt),
            "-A", "hydrogens",   # add hydrogens if missing
        ]

        print(f"Running MGLTools ligand preparation: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"MGLTools prepare_ligand4.py failed (exit {result.returncode}):\n"
                f"  stdout: {result.stdout.strip()}\n"
                f"  stderr: {result.stderr.strip()}"
            )

        if not out_pdbqt.exists():
            raise FileNotFoundError(
                f"MGLTools did not produce the expected output file: {out_pdbqt}\n"
                f"  stdout: {result.stdout.strip()}\n"
                f"  stderr: {result.stderr.strip()}"
            )

        pdbqt_string = out_pdbqt.read_text()

    return pdbqt_string


def main():
    parser = argparse.ArgumentParser(
        description="Unified script for Blind and Targeted molecular docking using Vina and MGLTools."
    )

    # Required Arguments
    parser.add_argument('-r', '--receptor', required=True, help="Path to the receptor PDBQT file")
    parser.add_argument('-l', '--ligand', required=True, help="Path to the ligand SDF file")
    parser.add_argument('-o', '--out_dir', required=True, help="Directory to save the resulting PDBQT poses")
    parser.add_argument('--prepare_receptor_script', default='/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py', help="Path to prepare_receptor4.py")

    # MGLTools Arguments
    parser.add_argument(
        '--mgltools_python',
        default='/opt/mgltools/bin/pythonsh',
        help="Path to the MGLTools Python 2 interpreter (pythonsh). Default: /opt/mgltools/bin/pythonsh"
    )
    parser.add_argument(
        '--prepare_ligand_script',
        default='/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py',
        help="Path to MGLTools' prepare_ligand4.py. Default: /opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    )

    # Optional Arguments for Targeted Docking
    parser.add_argument('--center', type=float, nargs=3,
                        help="Center coordinates (X Y Z). If omitted, blind docking is triggered.")
    parser.add_argument('--box_size', type=float, nargs=3, default=[15.0, 15.0, 15.0],
                        help="Box dimensions (X Y Z). Default: 15 15 15. Used only if --center is provided.")

    # Optional Vina Parameters
    parser.add_argument('--exhaustiveness', type=int, default=64,
                        help="Exhaustiveness of the search (default: 64)")
    parser.add_argument('--n_poses', type=int, default=20,
                        help="Number of poses to output (default: 20)")

    args = parser.parse_args()

    # 1. Setup paths and output directory
    os.makedirs(args.out_dir, exist_ok=True)
    rec_path = Path(args.receptor)
    lig_path = Path(args.ligand)

    # Extract base names for dynamic file naming (e.g., '6T0Y' and 'BADGE')
    rec_name = rec_path.stem.replace('_clean_h_meeko', '')
    lig_name = lig_path.stem

    # 2. Determine Mode & Grid Box Parameters
    if args.center:
        mode = "targeted"
        center_coords = args.center
        box_dims = args.box_size
        coords_str = f"{center_coords[0]:.1f}_{center_coords[1]:.1f}_{center_coords[2]:.1f}"
        out_filename = f"{rec_name}_{lig_name}_{mode}_{coords_str}.pdbqt"
    else:
        mode = "blind"
        print(f"Calculating grid box for {mode} docking...")
        center_coords, box_dims = get_receptor_center_and_size(args.receptor)
        out_filename = f"{rec_name}_{lig_name}_{mode}.pdbqt"

    out_filepath = os.path.join(args.out_dir, out_filename)

    print("-" * 40)
    print(f"Mode:           {mode.capitalize()} Docking")
    print(f"Receptor:       {args.receptor}")
    print(f"Ligand:         {args.ligand}")
    print(f"Center Coords:  {[round(c, 3) for c in center_coords]}")
    print(f"Box Dimensions: {[round(b, 3) for b in box_dims]}")
    print(f"Output File:    {out_filepath}")
    print("-" * 40)

    # 3. Initialize Vina & Receptor
    v = Vina(sf_name='vina')

    # v.set_receptor(str(rec_path))
    print("Preparing receptor with MGLTools...")
    prepared_receptor = prepare_receptor_mgltools(
        input_path=rec_path,
        mgltools_python=args.mgltools_python,
        prepare_receptor_script=args.prepare_receptor_script,
    )
    v.set_receptor(str(prepared_receptor))
    

    # 4. Prepare Ligand (MGLTools via subprocess)
    print("Preparing ligand with MGLTools (Python 2 subprocess)...")
    pdbqt_string = prepare_ligand_mgltools(
        sdf_path=lig_path,
        mgltools_python=args.mgltools_python,
        prepare_ligand_script=args.prepare_ligand_script,
    )

    # 5. Compute Maps & Set Ligand
    print("Computing Vina maps...")
    v.compute_vina_maps(center=center_coords, box_size=box_dims)
    v.set_ligand_from_string(pdbqt_string)

    # 6. Run Docking
    print(f"Docking (Exhaustiveness: {args.exhaustiveness}, Poses: {args.n_poses})...")
    v.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)

    # 7. Write Results
    v.write_poses(out_filepath, n_poses=args.n_poses, overwrite=True)
    print(f"Docking complete. Results saved to {out_filepath}")


if __name__ == "__main__":
    main()