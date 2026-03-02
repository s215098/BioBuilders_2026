import streamlit as st
import subprocess

### Functions ###

def msa_mafft(fasta_file):
    '''Perform multiple sequence alignment using MAFFT.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences to be aligned.
    Returns:
        fasta_file_aligned.fasta
    '''
    
    # use subprocess to call MAFFT in python
    
    result = subprocess.run(['mafft', '--auto', fasta_file], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed with error: {result.stderr}")

    # return the aligned sequences as a string which can be downloaded by the user
    return result.stdout

def build_tree(alignment_file):
    
    # convert to newick format using seqconverter
    # seqconverter --informat [input format] --outformat [output format] -i [alignment file name] > [output format file name].
    
    #
    return "tree.nwk"




### UI ###
st.title("Phylogenetic Clustering")

st.subheader("Multiple Sequence Alignment")
st.markdown("### Upload your sequences in FASTA format to perform multiple sequence alignment and build a phylogenetic tree.")
input_fasta = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "faa", "fna"])

if input_fasta is not None:
    st.success("FASTA file uploaded successfully!")
    st.write("File name:", input_fasta.name)
    st.write("File size:", input_fasta.size, "bytes")

    if st.button("Perform MSA"):
        # Call the function to read the FASTA file and perform MSA
        mafft_output = msa_mafft(input_fasta)
        st.success("Multiple sequence alignment completed successfully!")
        st.download_button("Download Alignment", data=mafft_output, file_name="alignment.fasta", mime="text/plain")
        # alignment = msa_mafft(input_fasta)

        # st.success("Multiple sequence alignment completed successfully!")
        # st.download_button("Download Alignment", data=alignment, file_name="alignment.fasta", mime="text/plain")

st.subheader("Phylogenetic Tree Construction")
st.markdown("### After performing multiple sequence alignment, you can build a phylogenetic tree to visualize the evolutionary relationships between the sequences.")
st.button("Build Phylogenetic Tree")

if st.button("Build Phylogenetic Tree"):
    # Call the functions to read the FASTA file, perform MSA, and build the tree
    fasta_file = input_fasta  # This should be the path to the uploaded FASTA file
    # sequences = read_fasta(fasta_file)
    # alignment = msa_mafft(fasta_file)
    # tree_file = build_tree(alignment)

    # st.success("Phylogenetic tree built successfully!")
    # st.download_button("Download Tree", data=open(tree_file, 'r').read(), file_name="phylo_tree.nwk", mime="text/plain")