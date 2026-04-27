#!/opt/anaconda3/envs/msa_env/bin/python

import streamlit as st

pages = {
    "About Us": [
        st.Page("biobuilder_info_app/biobuilder_info.py", title = "BioBuilder info", icon = "👩‍🔬"),
        st.Page("our_project_app/our_project.py", title = "Our project", icon = "🔬"),
    ],
    "Software": [
        st.Page("enzyme_search_app/enzyme_search.py", title = "Enzyme search", icon = "🔍"),
        # st.Page("substrate_search_app/substrate_search.py", title = "Substrate search", icon = "🧪"),
        st.Page("docking_app/docking.py", title = "Docking", icon = "⚓"),
        st.Page("phylo_clustering_app/phylo_clustering.py", title = "Phylogenetic clustering", icon = "🌳"),
    ],
}

pg = st.navigation(pages)
pg.run()
