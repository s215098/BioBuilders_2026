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
        
    log("Running Boltz-2 AI Prediction (Fetching MSA & Affinity)...")
    
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
        
    # Read and Log the AI Affinity
    if affinity_files:
        with open(affinity_files[0], 'r') as f:
            aff_data = json.load(f)
            log10_ic50 = aff_data.get("affinity_pred_value", 0.0)
            log(f"Boltz Predicted Affinity [log10(IC50)]: {log10_ic50:.2f}")
            
    return pdb_path, confidence

def main():
    parser = argparse.ArgumentParser(description="Single-Target Epoxy-Degradation Pipeline")
    parser.add_argument('-n', '--name', required=True, help="Name of the enzyme")
    parser.add_argument('-seq', '--sequence', required=True, help="Amino acid sequence")
    parser.add_argument('-smi', '--smiles', required=True, help="SMILES string")
    
    args = parser.parse_args()
    out_dir = f"outputs_boltz/{args.name}"
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"\n{'='*50}\n Target: {args.name}\n{'='*50}")
    
    try:
        boltz_pdb, ai_score = run_boltz(args.sequence, args.smiles, out_dir)
        log(f"Boltz Confidence: {ai_score:.3f}")
        
    except Exception as e:
        log(f"Pipeline Error: {str(e)}")

if __name__ == "__main__":
    main()
