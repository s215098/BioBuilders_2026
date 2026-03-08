"""
phylo_app.py
============
Streamlit app for:
  1. Multiple Sequence Alignment (MSA) via MAFFT
  2. Phylogenetic tree inference via FastTree  ← replaces PhyML (no binary install needed)
  3. Tree visualisation via Bio.Phylo

WHY FastTree instead of PhyML?
  - PhyML requires a system binary that is often not on PATH in cloud/server environments.
  - FastTree is available as a Python package (`pip install fasttree`) AND as a binary,
    but we call it via the `fasttree` Python wrapper so no manual PATH setup is needed.
  - FastTree still computes a proper maximum-likelihood tree (WAG/GTR model).

Install requirements (once, in your environment):
    pip install streamlit biopython fasttree matplotlib
    conda install -c bioconda mafft      # or: brew install mafft / apt install mafft
"""

# ── Standard library ──────────────────────────────────────────────────────────
import subprocess
import tempfile
import os
from io import StringIO
from typing import Optional, Tuple

# ── Third-party ───────────────────────────────────────────────────────────────
import streamlit as st

# Biopython: alignment I/O and tree drawing
try:
    from Bio import AlignIO, Phylo
    from Bio.Seq import Seq
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

# Detect FastTree binary on PATH (installed via brew/conda, not pip)
from shutil import which
FASTTREE_BIN = which("FastTree") or which("fasttree") or which("FastTreeMP")
FASTTREE_OK = FASTTREE_BIN is not None

# matplotlib: used by Bio.Phylo.draw() to render the tree
try:
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend (required for Streamlit)
    import matplotlib.pyplot as plt
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False



# ══════════════════════════════════════════════════════════════════════════════
#  HELPER: binary guard
# ══════════════════════════════════════════════════════════════════════════════

def _require_binary(name: str):
    """
    Check that `name` exists on PATH.
    Raises EnvironmentError with a helpful install hint if missing.
    """
    from shutil import which
    if which(name) is None:
        raise EnvironmentError(
            f"Required binary '{name}' was not found on your PATH.\n"
            f"  • macOS:  brew install {name}\n"
            f"  • Linux:  sudo apt install {name}  OR  conda install -c bioconda {name}\n"
            f"Then restart the app."
        )



# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1 – Multiple Sequence Alignment with MAFFT
# ══════════════════════════════════════════════════════════════════════════════

def assert_uniform_lengths(aligned_fasta: str):
    """
    Validate that every sequence in an *aligned* (gapped) FASTA has the
    same length.  Raises ValueError if lengths differ.
    """
    lengths = []
    seq_parts = []
    header_seen = False

    for line in aligned_fasta.splitlines():
        if line.startswith(">"):
            if header_seen:
                # Save the finished sequence length
                lengths.append(len("".join(seq_parts)))
            header_seen = True
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())

    # Don't forget the last sequence
    if header_seen:
        lengths.append(len("".join(seq_parts)))

    unique = set(lengths)
    if len(unique) != 1:
        raise ValueError(
            f"Alignment sanity-check failed – sequences have different lengths: {sorted(unique)}"
        )


def run_mafft(uploaded_file) -> str:
    """
    Align sequences in `uploaded_file` (a Streamlit UploadedFile object)
    using MAFFT --auto.

    Returns:
        Aligned FASTA as a plain string.

    Raises:
        EnvironmentError  – if `mafft` binary is missing
        RuntimeError      – if MAFFT exits with a non-zero return code
        ValueError        – if the resulting alignment has inconsistent lengths
    """
    _require_binary("mafft")           # die early with a clear message

    uploaded_file.seek(0)              # reset file pointer (may have been read before)

    # Write the uploaded bytes to a temporary file on disk so MAFFT can read it
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["mafft", "--auto", tmp_path],
            capture_output=True,
            text=True
        )
    finally:
        os.unlink(tmp_path)            # clean up temp file even if subprocess fails

    if result.returncode != 0:
        raise RuntimeError(f"MAFFT exited with code {result.returncode}:\n{result.stderr}")

    aligned = result.stdout

    # Tidy up: remove blank lines, normalise gap character
    aligned = "\n".join(line for line in aligned.splitlines() if line.strip())
    aligned = aligned.replace(".", "-")   # some aligners use '.' for gaps

    assert_uniform_lengths(aligned)        # sanity check before returning
    return aligned


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2a – Convert aligned FASTA → PHYLIP (relaxed) via Biopython
# ══════════════════════════════════════════════════════════════════════════════

def fasta_to_phylip(aligned_fasta: str) -> str:
    """
    Convert an aligned FASTA string to PHYLIP-relaxed format using Biopython.

    PHYLIP-relaxed keeps full sequence names (classic PHYLIP truncates to 10
    characters, which is almost always a problem for real data).

    Returns:
        PHYLIP text as a plain string.
    """
    if not BIOPYTHON_OK:
        raise ImportError("Biopython is not installed. Run: pip install biopython")

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_path  = os.path.join(tmpdir, "aligned.fasta")
        phy_path = os.path.join(tmpdir, "aligned.phy")

        # Write FASTA to disk so AlignIO can parse it
        with open(fa_path, "w") as fh:
            fh.write(aligned_fasta)

        aln = AlignIO.read(fa_path, "fasta")

        # Normalise each sequence: uppercase, gap character '-'
        for rec in aln:
            rec.seq = Seq(str(rec.seq).upper().replace(".", "-").replace(" ", ""))

        # Confirm lengths are still uniform after normalisation
        lengths = {len(rec.seq) for rec in aln}
        if len(lengths) != 1:
            raise ValueError(f"Alignment has unequal lengths after normalisation: {lengths}")

        AlignIO.write(aln, phy_path, "phylip-relaxed")

        with open(phy_path) as fh:
            phylip_text = fh.read()

    if not phylip_text.strip():
        raise RuntimeError("PHYLIP conversion produced empty output – check your alignment.")

    return phylip_text


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2b – Build ML tree with FastTree
# ══════════════════════════════════════════════════════════════════════════════

def run_fasttree(aligned_fasta: str, is_protein: bool = True) -> str:
    """
    Run FastTree on an aligned FASTA string and return a Newick tree.

    FastTree is a maximum-likelihood tree builder that is:
      • Much faster than PhyML for large alignments
      • Available as a Python package (pip install fasttree) – no PATH fiddling
      • Defaults to WAG+CAT model for proteins, GTR+CAT for nucleotides

    Args:
        aligned_fasta:  Aligned sequences in FASTA format (string).
        is_protein:     True → amino-acid sequences; False → nucleotide sequences.

    Returns:
        Newick tree string.
    """
    
    if not FASTTREE_OK:
        raise EnvironmentError(
            "FastTree not found. Install with: brew install fasttree"
        )
    cmd = [FASTTREE_BIN, "-quiet"]

    if is_protein:
        cmd += ["-wag"]
    else:
        cmd += ["-nt", "-gtr"]
    result = subprocess.run(cmd, input=aligned_fasta, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"FastTree failed:\n{result.stderr}")
    newick = result.stdout.strip()

    if not newick:
        raise RuntimeError(f"FastTree produced no output.\n{result.stderr}")
    
    return newick


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3 – Visualise Newick tree with Bio.Phylo + matplotlib
# ══════════════════════════════════════════════════════════════════════════════

def draw_tree(newick: str):
    """
    Render a Newick tree as a polished matplotlib figure.

    Improvements over the default Bio.Phylo.draw():
      - Warm off-white background instead of near-black (much easier to read)
      - Thicker, dark teal branch lines with high contrast
      - Larger fonts throughout (labels, axes, title)
      - Confidence values hidden (they overlap taxon names and clutter the plot)
      - Clean spines and subtle grid on the x-axis only
      - Extra left margin so long taxon names don't get clipped
    """
    if not BIOPYTHON_OK or not MATPLOTLIB_OK:
        return None

    tree = Phylo.read(StringIO(newick), "newick")
    tree.ladderize()

    n_leaves   = tree.count_terminals()
    fig_height = max(6, min(n_leaves * 0.55, 40))   # taller per leaf than before

    BG      = "#f7f4ef"   # warm off-white background
    BRANCH  = "#1a6b72"   # dark teal — high contrast on the warm bg
    TEXT    = "#1c2b2d"   # near-black text
    SPINE   = "#c5bfb0"   # soft warm grey for axes

    fig, ax = plt.subplots(figsize=(13, fig_height))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    # Bio.Phylo.draw places confidence values as text objects mixed with taxon
    # labels.  We suppress them by temporarily zeroing out confidence values,
    # then restoring after drawing — cleaner than post-hoc text deletion.
    for clade in tree.find_clades():
        clade._orig_confidence = clade.confidence
        clade.confidence = None     # hide confidence values on the plot

    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    # Restore confidence values (in case the tree object is reused)
    for clade in tree.find_clades():
        clade.confidence = clade._orig_confidence

    # ── Restyle branch lines ──────────────────────────────────────────────────
    for line in ax.get_lines():
        line.set_color(BRANCH)
        line.set_linewidth(1.8)

    # ── Restyle all text (taxon labels + any remaining annotations) ───────────
    for text_obj in ax.texts:
        text_obj.set_color(TEXT)
        text_obj.set_fontsize(12)
        text_obj.set_fontfamily("monospace")

    # ── Axes styling ──────────────────────────────────────────────────────────
    ax.set_xlabel("branch length", fontsize=12, color=TEXT, labelpad=8)
    ax.set_ylabel("taxa", fontsize=12, color=TEXT, labelpad=8)
    ax.tick_params(axis="x", colors=TEXT, labelsize=10)
    ax.tick_params(axis="y", colors=TEXT, labelsize=10)

    # Subtle x-axis grid only (branch length reference lines)
    ax.xaxis.grid(True, color=SPINE, linewidth=0.5, linestyle="--", alpha=0.6)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_color(SPINE)
        spine.set_linewidth(0.8)

    # Extra left margin so long names aren't clipped
    ax.margins(x=0.02)
    fig.subplots_adjust(left=0.22)

    ax.set_title(
        "Phylogenetic Tree (FastTree ML)",
        fontsize=15, fontweight="bold",
        color=TEXT, pad=14
    )

    plt.tight_layout()
    return fig



# ══════════════════════════════════════════════════════════════════════════════
#  STREAMLIT UI
# ══════════════════════════════════════════════════════════════════════════════

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Phylo Suite",
    page_icon="🌳",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Custom CSS for a clean dark-lab aesthetic ─────────────────────────────────
st.markdown("""
<style>
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

  html, body, [class*="css"] {
      font-family: 'IBM Plex Sans', sans-serif;
      background-color: #0f1923;
      color: #cfd8dc;
  }
  h1, h2, h3 { font-family: 'IBM Plex Mono', monospace; color: #4fc3f7; letter-spacing: -0.5px; }
  .stButton > button {
      background: #0d47a1;
      color: #e3f2fd;
      border: none;
      border-radius: 4px;
      font-family: 'IBM Plex Mono', monospace;
      font-weight: 600;
      padding: 0.5rem 1.4rem;
      transition: background 0.2s;
  }
  .stButton > button:hover { background: #1565c0; }
  .stDownloadButton > button {
      background: #004d40;
      color: #b2dfdb;
      border: none;
      border-radius: 4px;
      font-family: 'IBM Plex Mono', monospace;
  }
  .stDownloadButton > button:hover { background: #00695c; }
  .stAlert { border-radius: 6px; }
  code { font-family: 'IBM Plex Mono', monospace; color: #80cbc4; }
  hr { border-color: #1e3a4a; }
  .step-badge {
      display: inline-block;
      background: #0d47a1;
      color: #e3f2fd;
      font-family: 'IBM Plex Mono', monospace;
      font-size: 0.75rem;
      padding: 2px 10px;
      border-radius: 12px;
      margin-bottom: 8px;
  }
</style>
""", unsafe_allow_html=True)

# ── Title ─────────────────────────────────────────────────────────────────────
st.markdown("# 🌳 Phylo Suite")
st.markdown("This app carries out Multiple Sequence Alignment & produces a maximum-likelihood phylogenetic tree in two clicks.")
st.markdown("---")


# ── Dependency check banner ───────────────────────────────────────────────────
missing = []
if not BIOPYTHON_OK:
    missing.append("`pip install biopython`")
if not FASTTREE_OK:
    missing.append("`pip install fasttree`")
if not MATPLOTLIB_OK:
    missing.append("`pip install matplotlib`")

if missing:
    st.warning(
        "⚠️ Some dependencies are missing. Install them and restart:\n\n"
        + "\n".join(missing)
    )


# ═════════════════════════════════════════
#  STEP 1 – Upload + MSA
# ═════════════════════════════════════════
st.markdown('<div class="step-badge">STEP 1</div>', unsafe_allow_html=True)
st.subheader("Multiple Sequence Alignment")
st.markdown("Upload a FASTA file (protein or nucleotide). MAFFT will align the sequences.")

col_upload, col_info = st.columns([3, 1])

with col_upload:
    input_fasta = st.file_uploader(
        "Upload FASTA file",
        type=["fasta", "fa", "faa", "fna"],
        help="Accepted extensions: .fasta, .fa, .faa, .fna"
    )

with col_info:
    if input_fasta:
        st.success(f"✓ {input_fasta.name}")
        st.caption(f"Size: {input_fasta.size:,} bytes")

# ── Run MSA ───────────────────────────────────────────────────────────────────
if input_fasta:
    if st.button("▶ Run MAFFT Alignment", key="run_msa"):
        with st.spinner("Aligning sequences with MAFFT…"):
            try:
                aligned = run_mafft(input_fasta)
                st.session_state["alignment"] = aligned   # store for the next step
                st.success("✓ Alignment complete!")
            except Exception as e:
                st.error(f"MAFFT error:\n\n{e}")
                st.session_state.pop("alignment", None)   # clear stale result

    # vertical space
    st.markdown("")

    # Show download button only if we have a valid alignment in session
    if "alignment" in st.session_state:
        st.download_button(
            "⬇ Download Alignment (.fasta)",
            data=st.session_state["alignment"].encode(),
            file_name="alignment.fasta",
            mime="text/plain",
        )

        # vertical space
        st.markdown("")

        # Show a small preview (first 500 chars) so the user can sanity-check
        with st.expander("Preview aligned FASTA"):
            st.code(st.session_state["alignment"][:800] + " …", language="text")


# ═════════════════════════════════════════
#  STEP 2 – Build phylogenetic tree
# ═════════════════════════════════════════
st.markdown("---")
st.markdown('<div class="step-badge">STEP 2</div>', unsafe_allow_html=True)
st.subheader("Phylogenetic Tree Construction")

# ── Options ───────────────────────────────────────────────────────────────────
col_phyl_main, col_phyl_info = st.columns([1.5, 1])

with col_phyl_main:
    st.markdown(
    "FastTree builds a maximum-likelihood tree from the alignment above.\n"
    )
    seq_type = st.radio(
        "Sequence type",
        options=["Protein (amino acid)", "Nucleotide (DNA/RNA)"],
        index=0,
        help="""Protein → Whelan and Goldman model.  
                Nucleotide → General Time Reversible model.""",
    )
    is_protein = seq_type.startswith("Protein")

with col_phyl_info:
    st.info(
        "**About the model**\n\n"
        f"{'Protein: Whelan and Goldman + Category substitution model' if is_protein else 'Nucleotides: General Time Reversible + Category substitution model'}\n\n"
        """\nThe program uses "FastTree" which utilizes the fast approximate-ML and is well-suited for both small and large datasets."""
    )

# ── Run tree ──────────────────────────────────────────────────────────────────
if st.button("▶ Build Phylogenetic Tree", key="run_tree"):
    if "alignment" not in st.session_state:
        st.error("Please complete Step 1 (MSA) first.")
    else:
        with st.spinner("Running FastTree… (may take a minute for large datasets)"):
            try:
                newick = run_fasttree(
                    st.session_state["alignment"],
                    is_protein=is_protein
                )
                st.session_state["newick"] = newick
                st.success("✓ Tree built successfully!")
            except Exception as e:
                st.error(f"Tree-building error:\n\n{e}")
                st.session_state.pop("newick", None)

# ── Display results ───────────────────────────────────────────────────────────
if "newick" in st.session_state:
    newick = st.session_state["newick"]

    # vertical space
    st.markdown("")
    
    # Download Newick
    st.download_button(
        "⬇ Download Newick Tree (.nwk)",
        data=newick.encode(),
        file_name="tree_fasttree.nwk",
        mime="text/plain",
    )

    # vertical space
    st.markdown("")
    
    # Show raw Newick text
    with st.expander("Show raw Newick string"):
        st.code(newick, language="text")

    # vertical space
    st.markdown("")
    
    # Draw tree graphic
    if MATPLOTLIB_OK and BIOPYTHON_OK:
        with st.spinner("Rendering tree…"):
            try:
                fig = draw_tree(newick)
                if fig:
                    st.pyplot(fig)
                    plt.close(fig)     # free memory
            except Exception as e:
                st.warning(f"Could not render tree graphic: {e}\nYou can still download the Newick file above.")
    else:
        st.info("Install matplotlib and biopython to see a tree graphic here.")

# ── Footer ────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption(
    "Phylo Suite · MSA via [MAFFT](https://mafft.cbrc.jp/alignment/software/) · "
    "Tree via [FastTree](http://www.microbesonline.org/fasttree/) · "
    "Built with [Streamlit](https://streamlit.io) & [Biopython](https://biopython.org)"
)