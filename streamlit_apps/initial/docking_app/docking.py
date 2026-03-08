import streamlit as st
import tempfile
import os
import subprocess
import numpy as np
import py3Dmol
from stmol import showmol
from vina import Vina
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# --- HELPER FUNCTIONS ---

def get_receptor_center_and_size(pdbqt_file):
    """Calculates bounding box from a PDB/PDBQT file."""
    coordinates = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinates.append([x, y, z])
    
    coords = np.array(coordinates)
    center = np.mean(coords, axis=0) 
    size = np.max(coords, axis=0) - np.min(coords, axis=0)
    size = size + 10.0 
    
    return center.tolist(), size.tolist()

def parse_vina_poses(pdbqt_str):
    """Splits a multi-model PDBQT string and extracts affinity scores."""
    poses = []
    scores = []
    current_pose = []
    for line in pdbqt_str.splitlines():
        current_pose.append(line)
        
        if line.startswith("REMARK VINA RESULT:"):
            # Extract the affinity score
            score = float(line.split()[3])
            scores.append(score)
            
        if line.strip() == "ENDMDL":
            poses.append("\n".join(current_pose))
            current_pose = []
            
    return poses, scores

# --- STREAMLIT APP ---

st.title("Molecular Docking Pipeline")

# 1. Inputs
st.subheader("1. Upload Inputs")
col1, col2 = st.columns(2)
with col1:
    receptor_file = st.file_uploader("Upload Receptor (PDB)", type=["pdb"])
with col2:
    ligand_file = st.file_uploader("Upload Ligand (SDF)", type=["sdf"])

# Exposing your hardcoded parameters to the UI for flexibility
st.subheader("2. Settings")
col3, col4 = st.columns(2)
with col3:
    exhaustiveness = st.slider("Exhaustiveness", min_value=1, max_value=128, value=64)
with col4:
    n_poses = st.number_input("Number of Poses", min_value=1, max_value=20, value=20)

# 2. Execution Block
if st.button("Run Docking Pipeline") and receptor_file and ligand_file:
    with st.spinner("Preparing files and running Vina..."):
        
        with tempfile.TemporaryDirectory() as tmpdir:
            raw_pdb = os.path.join(tmpdir, "receptor.pdb")
            rec_basename = os.path.join(tmpdir, "receptor_meeko")
            rec_pdbqt = f"{rec_basename}.pdbqt"
            
            lig_sdf = os.path.join(tmpdir, "ligand.sdf")
            out_pdbqt = os.path.join(tmpdir, "docking_results.pdbqt")

            # Write uploaded files to disk
            with open(raw_pdb, "wb") as f: f.write(receptor_file.getvalue())
            with open(lig_sdf, "wb") as f: f.write(ligand_file.getvalue())

            try:
                # Step A: Receptor Prep (Using your exact requested command)
                st.text("Converting receptor via mk_prepare_receptor.py...")
                subprocess.run([
                    "mk_prepare_receptor.py", 
                    "-i", raw_pdb, 
                    "-o", rec_basename, 
                    "-p"
                ], check=True, capture_output=True, text=True)

                # Save the raw PDB text for visualization later
                with open(raw_pdb, "r") as f:
                    st.session_state['receptor_pdb_str'] = f.read()

                # Step B: Box Dimensions
                center_coords, box_dims = get_receptor_center_and_size(rec_pdbqt)

                # Step C: Ligand Prep
                st.text("Optimizing ligand 3D structure with RDKit...")
                suppl = Chem.SDMolSupplier(lig_sdf)
                mol = next(suppl)
                mol_with_h = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol_with_h)
                
                preparer = MoleculePreparation()
                preparer.prepare(mol_with_h)
                ligand_pdbqt_string = preparer.write_pdbqt_string()

                # Step D: AutoDock Vina Execution
                st.text("Computing Vina maps and running docking...")
                v = Vina(sf_name='vina')
                v.set_receptor(rec_pdbqt)
                v.compute_vina_maps(center=center_coords, box_size=box_dims)
                v.set_ligand_from_string(ligand_pdbqt_string)
                
                v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
                v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)

                # Read final poses into memory
                with open(out_pdbqt, "r") as f:
                    st.session_state['docked_poses_str'] = f.read()

                st.success("Docking complete!")

            except subprocess.CalledProcessError as e:
                st.error("Meeko receptor preparation failed.")
                st.code(e.stderr)
            except Exception as e:
                st.error(f"An error occurred: {e}")

# --- VISUALIZATION BLOCK ---
if 'docked_poses_str' in st.session_state and 'receptor_pdb_str' in st.session_state:
    st.subheader("3. Interactive Results Viewer")
    
    # Parse the poses and scores from the Vina output
    poses, scores = parse_vina_poses(st.session_state['docked_poses_str'])
    
    col_viz1, col_viz2 = st.columns([1, 3]) 
    
    with col_viz1:
        st.write("### Controls")
        # Create a dropdown menu showing the pose number and its affinity score
        pose_options = [f"Pose {i+1} ({scores[i]} kcal/mol)" for i in range(len(scores))]
        selected_pose_label = st.selectbox("Select Docking Pose:", pose_options)
        
        # Determine which pose index the user selected
        pose_idx = pose_options.index(selected_pose_label)
        
        # Display the affinity prominently
        st.metric(label="Binding Affinity", value=f"{scores[pose_idx]} kcal/mol")
        
        st.download_button(
            label="Download All Poses (PDBQT)",
            data=st.session_state['docked_poses_str'],
            file_name="vina_results.pdbqt",
            mime="chemical/x-pdbqt"
        )
        
    with col_viz2:
        # Render the specific pose chosen by the user
        view = py3Dmol.view(width=600, height=500)
        
        # Add receptor (using the original PDB so it looks nice)
        view.addModel(st.session_state['receptor_pdb_str'], 'pdb')
        view.setStyle({'model': 0}, {'cartoon': {'color': 'white'}})
        
        # Add ONLY the selected ligand pose (using the raw PDBQT)
        view.addModel(poses[pose_idx], 'pdbqt')
        view.setStyle({'model': 1}, {
            'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}, 
            'sphere': {'colorscheme': 'cyanCarbon', 'radius': 0.4}
        })
        
        view.zoomTo({'model': 1}) 
        
        showmol(view, height=500, width=600)