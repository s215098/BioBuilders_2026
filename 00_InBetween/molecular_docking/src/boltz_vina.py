import os
import sys
import json
import subprocess
import argparse


def log(msg):
    print(f"[EPOXY-SCREEN] {msg}")

def run_boltz(sequence, smiles, out_dir):
    log("Generating Boltz YAML input file...")
    yaml_path = f"{out_dir}/boltz_input.yaml"
    
    yaml_content = f"""
sequences:
  - protein:
      id: A
      sequence: {sequence}
  - ligand:
      id: B
      smiles: '{smiles}'
properties:
  - affinity:
      binder: B
"""
    with open(yaml_path, "w") as f:
        f.write(yaml_content.strip())
        
    log("Running Boltz-2 Prediction (Fetching MSA & Affinity)...")
    
    subprocess.run([
        "boltz", "predict", 
        yaml_path, 
        "--out_dir", f"{out_dir}/boltz_out",
        "--use_msa_server",
        "--override"
    ], check=True, stdout=subprocess.DEVNULL) 
    
    # File finder: Boltz outputs are nested and we need to find the confidence, affinity, and structure files.

    import glob
    boltz_out_dir = f"{out_dir}/boltz_out"
    
    # 1. Find Confidence JSON
    json_files = glob.glob(f"{boltz_out_dir}/**/confidence_*.json", recursive=True)
    if not json_files:
        raise FileNotFoundError(f"Boltz finished, but could not find confidence JSON in {boltz_out_dir}")
    json_path = json_files[0]
    
    # 2. Find Affinity JSON
    affinity_files = glob.glob(f"{boltz_out_dir}/**/affinity_*.json", recursive=True)
    
    # 3. Find Structure (.cif or .pdb)
    cif_files = glob.glob(f"{boltz_out_dir}/**/*.cif", recursive=True)
    pdb_files = glob.glob(f"{boltz_out_dir}/**/*.pdb", recursive=True)
    
    if cif_files:
        pdb_path = cif_files[0]
    elif pdb_files:
        pdb_path = pdb_files[0]
    else:
        raise FileNotFoundError(f"Boltz finished, but could not find structure file in {boltz_out_dir}")
    
    # Read the Confidence
    with open(json_path, 'r') as f:
        data = json.load(f)
        confidence = data.get("confidence_score", data.get("ptm", 0.0))
        
    # Read and Log the Affinity
    if affinity_files:
        with open(affinity_files[0], 'r') as f:
            aff_data = json.load(f)
            log10_ic50 = aff_data.get("affinity_pred_value", 0.0)
            log(f"Boltz Predicted Affinity [log10(IC50)]: {log10_ic50:.2f}")
            
    return pdb_path, confidence

def extract_pocket_and_clean(pdb_file, clean_out_path):
    import numpy as np
    from Bio.PDB import PDBParser, PDBIO, Select

    class ReceptorOnlySelect(Select):
        def accept_residue(self, residue):
            if residue.get_resname() in ["UNL", "LIG1", "HOH"]:
                return 0 
            return 1

    log("Extracting active site and stripping AI ligand...")
    parser = PDBParser(QUIET=True)
    
    if pdb_file.endswith(".cif"):
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
        
    structure = parser.get_structure("complex", pdb_file)
    coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in ["LIG1"]:
                    for atom in residue:
                        coords.append(atom.get_coord())
                        
    if not coords:
        raise ValueError("Ligand LIG1 not found in Boltz prediction.")
        
    coords = np.array(coords)
    center = coords.mean(axis=0).tolist()
    
    raw_span = coords.max(axis=0) - coords.min(axis=0)
    box_size = np.maximum(raw_span + 10.0, [20.0, 20.0, 20.0]).tolist()
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_out_path, ReceptorOnlySelect())
    
    return center, box_size

def prep_receptor(clean_pdb, pdbqt_out):
    import subprocess
    import os
    
    log("Preparing Receptor PDBQT...")
    
    cmd = [
        "mk_prepare_receptor.py", 
        "--read_pdb", clean_pdb, 
        "-p", pdbqt_out
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError("Failed to prepare receptor.")
        
    if not os.path.exists(pdbqt_out):
        raise FileNotFoundError(f"Meeko ran, but failed to create {pdbqt_out}")
        
    log(f"Receptor PDBQT successfully saved to {pdbqt_out}")

def run_native_docking(receptor_pdbqt, ligand_sdf, center, box, out_file):
    from vina import Vina
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation

    log("Minimizing custom substrate with RDKit & Meeko...")
    suppl = Chem.SDMolSupplier(ligand_sdf)
    mol = next(suppl)
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_with_h)
    
    preparer = MoleculePreparation()
    preparer.prepare(mol_with_h)
    ligand_pdbqt_string = preparer.write_pdbqt_string()
    
    log(f"Executing Native Vina at coordinates {[round(c, 2) for c in center]}...")
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=center, box_size=box)
    v.set_ligand_from_string(ligand_pdbqt_string)
    
    v.dock(exhaustiveness=32, n_poses=5)
    v.write_poses(out_file, n_poses=5, overwrite=True)
    
    best_score = v.energies()[0][0] 
    return best_score

def main():
    parser = argparse.ArgumentParser(description="Single-Target Epoxy-Degradation Pipeline")
    parser.add_argument('-n', '--name', required=True, help="Name of the enzyme")
    parser.add_argument('-seq', '--sequence', required=True, help="Amino acid sequence")
    parser.add_argument('-sdf', '--substrate', required=True, help="SDF file")
    parser.add_argument('-smi', '--smiles', required=True, help="SMILES string")
    parser.add_argument('-t', '--threshold', type=float, default=0.60, help="Min Boltz confidence")
    
    args = parser.parse_args()
    out_dir = f"outputs/{args.name}"
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"\n{'='*50}\n Target: {args.name}\n{'='*50}")
    
    try:
        boltz_pdb, ai_score = run_boltz(args.sequence, args.smiles, out_dir)
        log(f"Boltz Confidence: {ai_score:.3f}")
        
        if ai_score < args.threshold:
            log(f"Failed AI Threshold (< {args.threshold}). Skipping physics docking.")
            sys.exit(0)
            
        clean_pdb = f"{out_dir}/receptor_clean.pdb"
        receptor_pdbqt = f"{out_dir}/receptor.pdbqt"
        center, box = extract_pocket_and_clean(boltz_pdb, clean_pdb)
        prep_receptor(clean_pdb, receptor_pdbqt)
        
        final_poses = f"{out_dir}/final_poses.pdbqt"
        vina_score = run_native_docking(receptor_pdbqt, args.substrate, center, box, final_poses)
        
        print(f"\n{'='*50}")
        print(f" PIPELINE SUCCESS ")
        print(f" Target: {args.name}")
        print(f" Vina Binding Affinity: {vina_score:.2f} kcal/mol")
        print(f" Best pose saved to: {final_poses}")
        print(f"{'='*50}\n")
        
    except Exception as e:
        log(f"Pipeline Error: {str(e)}")

if __name__ == "__main__":
    main()
