#!/opt/anaconda3/envs/msa_env/bin/python

import streamlit as st

pages = {
    "About Us": [
        st.Page("biobuilder_info.py", title="BioBuilder info"),
        st.Page("our_project.py", title="Our project"),
    ],
    "Software": [
        st.Page("enzyme_search.py", title="Enzyme search"),
        st.Page("substrate_search.py", title="Substrate search"),
    ],
}

pg = st.navigation(pages)
pg.run()
