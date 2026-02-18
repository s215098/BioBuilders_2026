import streamlit as st

st.title("Enzyme Search")

# some example content for the enzyme search page

# chose substrate from dropdown menu
substrate = st.selectbox("Select a substrate", ["Substrate A", "Substrate B", "Substrate C"])

# Example data
data = {
    "Substrate A": {"Enzyme 1": "High", "Enzyme 2": "Low", "Enzyme 3": "None"},
    "Substrate B": {"Enzyme 1": "Low", "Enzyme 2": "High", "Enzyme 3": "Low"},
    "Substrate C": {"Enzyme 1": "None", "Enzyme 2": "Low", "Enzyme 3": "High"},
}

# display results based on selected substrate
if substrate == "Substrate A":
    st.write("Enzyme 1 is known to interact with Substrate A.")
    # plot for this using example data
    st.bar_chart(data["Substrate A"])

    # Modelling plot example
    st.subheader("Modeling")
    st.image("streamlit_apps/initial/example_modelling_plot.png", caption="Modeling plot for Enzyme 1 and Substrate A")
elif substrate == "Substrate B":
    st.write("Enzyme 2 is known to interact with Substrate B.")
    st.bar_chart(data["Substrate B"])

    # Modelling plot example
    st.subheader("Modeling")
    st.image("streamlit_apps/initial/example_modelling_plot.png", caption="Modeling plot for Enzyme 2 and Substrate B")

elif substrate == "Substrate C":
    st.write("Enzyme 3 is known to interact with Substrate C.")
    st.bar_chart(data["Substrate C"])
    # Modelling plot example
    st.subheader("Modeling")
    st.image("streamlit_apps/initial/example_modelling_plot.png", caption="Modeling plot for Enzyme 3 and Substrate C")
