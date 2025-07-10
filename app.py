from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np

from plots import (
    plot_graph,
    plot_molecular_similarity_heatmap,
    plot_proteomics_barplots,
    plot_genomics_barplots,
    radar_plot,
)

# Constants
D_NAME = "DRUG_NAME"
G_DRIVER = "Driver_Gene"
G_NAME = "GENE"
D_ACT = "ACT_VALUE"

# Set page title and layout
st.set_page_config(page_title="Drug-Gene Interaction Explorer", layout="wide")
st.title("Drug-Gene Interaction Explorer")

# Load data
@st.cache_data
def load_data():
    PROCESSED_DATA_PATH = Path("data/DrugCentral/processed_data/")

    df_cancer_gene = pd.read_csv(PROCESSED_DATA_PATH / "gene.id.csv")
    df_gene_drug = pd.read_table(
        PROCESSED_DATA_PATH / "drug.target.interaction.fda.cosmic.cancer.type.tsv"
    )
    df_genomics = pd.read_csv(
        "data/LINCS/processed_data/lincs.drugcentral.csv"
    )
    df_proteomics = pd.read_csv(
        "data/DeepCoverageMOA/processed_data/deepcoveragemoa.lincs.csv"
    )

    # Load SMILES data for molecular similarity
    df_drug_to_smiles = pd.read_table(
        "data/DrugCentral/raw_data/structures.smiles.tsv"
    )[["INN", "SMILES"]]

    drug_name_to_smiles = {t.INN: t.SMILES for t in df_drug_to_smiles.itertuples()}

    return df_cancer_gene, df_gene_drug, df_genomics, df_proteomics, drug_name_to_smiles


# Load the data with caching
df_cancer_gene, df_gene_drug, df_genomics, df_proteomics, drug_name_to_smiles = load_data()

# Filter the gene-drug data to include only Kd and IC50 activity types
df_gene_drug = df_gene_drug.query("ACT_TYPE in ['Kd', 'IC50']")

# Sidebar for cancer selection
st.sidebar.header("Data Selection")
cancers_options = df_cancer_gene["CANCER"].unique()
cancer = st.sidebar.selectbox("Select cancer type:", cancers_options)

# Filter genes associated with the selected cancer
genes_in_cancer = df_cancer_gene.query(f"CANCER=='{cancer}'")["GENE"]
df_gene_drug_filtered = df_gene_drug.query("GENE in @genes_in_cancer")

# Show data overview
st.header(f"Cancer Type: {cancer}")
st.write(f"Number of genes associated with {cancer}: {genes_in_cancer.nunique()}")

# Display dataframe with custom styling
with st.expander("View Gene-Drug Interaction Data"):
    st.dataframe(
        df_gene_drug_filtered.style.highlight_max(axis=0, subset=[D_ACT]),
        use_container_width=True,
    )

# Visualizations
st.header("Visualizations")

# Create tabs for different visualizations
tab1, tab2 = st.tabs(["Drug Activity Radar", "Drug-Gene Network"])

with tab1:
    st.subheader("Drug Activity Radar Plot")
    st.write("This plot shows drug activity values grouped by activity type.")

    fig_radar = radar_plot(df_gene_drug_filtered)
    st.pyplot(fig_radar)

with tab2:
    st.subheader("Drug-Gene Network")
    st.write("This network visualization shows relationships between drugs and genes.")

    fig_network = plot_graph(df_gene_drug_filtered)
    st.pyplot(fig_network)

# Drug selection for similarity analysis
st.header("Drug Similarity Analysis")

# Filter to drugs that have genomics data available
drugs_with_genomics = set(df_genomics["DRUG_NAME_1"].unique()) & set(
    df_gene_drug_filtered[D_NAME].unique()
)

# Display number of drugs with genomics data
st.write(f"Number of drugs with genomics data: {len(drugs_with_genomics)}")

# Filter to drugs that have proteomics data available
drugs_with_proteomics = set(df_proteomics["DRUG_NAME_1"].unique()) & set(
    df_gene_drug_filtered[D_NAME].unique()
)

# Display number of drugs with proteomics data
st.write(f"Number of drugs with proteomics data: {len(drugs_with_proteomics)}")

# Filter to drugs that have genomics and proteomics data available
drugs_with_genomics_proteomics = set(list(drugs_with_genomics) + list(drugs_with_proteomics))

# Create a multi-select widget to select drugs with genomics data
selected_genomics_drugs = st.multiselect(
    "Select drugs to analyze genomic profile similarities:", sorted(drugs_with_genomics)
)

# If drugs are selected, show genomics and structural similarities
if selected_genomics_drugs:
    st.subheader("Genomics Similarity")
    st.write("These plots show the most similar drugs based on genomic profiles.")

    # Get Genomics similarity data and plots
    genomics_figs, similar_genomics_drugs = plot_genomics_barplots(df_genomics, selected_genomics_drugs)

    # Display genomics similarity plots
    for genomic_fig in genomics_figs:
        # st.pyplot(genomic_fig)
        st.plotly_chart(genomic_fig, use_container_width=True)

    # Display structural similarity if drugs have SMILES data
    st.subheader("Molecular Structure Similarity")
    st.write(
        "This heatmap shows Tanimoto similarity between drug molecular structures."
    )

    # Check if the selected drugs have SMILES data
    drugs_with_smiles = [d for d in selected_genomics_drugs if d in drug_name_to_smiles]
    similar_drugs_with_smiles = [d for d in similar_genomics_drugs if d in drug_name_to_smiles]

    if drugs_with_smiles and similar_drugs_with_smiles:
        fig = plot_molecular_similarity_heatmap(
            drugs_with_smiles, similar_drugs_with_smiles, drug_name_to_smiles
        )
        st.pyplot(fig)
    else:
        st.warning("No molecular structure data available for the selected drugs.")

    # Display raw similarity data in expander
    with st.expander("View similarity data"):
        # Create a dataframe of all pairwise similarities
        similarity_data = []

        for drug in selected_genomics_drugs:
            # Get similar drugs for this drug
            drug_data = df_genomics[df_genomics["DRUG_NAME_1"] == drug]

            if not drug_data.empty:
                # Sort by similarity score
                drug_data = drug_data.sort_values(
                    by="LINCS Pearson (r)", ascending=False
                )

                # Add to our collection
                similarity_data.append(drug_data)

        if similarity_data:
            combined_data = pd.concat(similarity_data)
            st.dataframe(combined_data, use_container_width=True)
        else:
            st.write("No similarity data available for the selected drugs.")

# Create a multi-select widget to select drugs with proteomics data
selected_proteomics_drugs = st.multiselect(
    "Select drugs to analyze proteomic profile similarities:", sorted(drugs_with_proteomics)
)

# If drugs are selected, show proteomics and structural similarities
if selected_proteomics_drugs:
    st.subheader("Proteomics Similarity")
    st.write("These plots show the most similar drugs based on proteomics profiles.")

    # Get proteomics similarity data and plots
    proteomics_figs, similar_proteomics_drugs = plot_proteomics_barplots(df_proteomics, selected_proteomics_drugs)

    # Display proteomics similarity plots
    for pro_fig in proteomics_figs:
        st.pyplot(pro_fig)

    # Display structural similarity if drugs have SMILES data
    st.subheader("Molecular Structure Similarity")
    st.write(
        "This heatmap shows Tanimoto similarity between drug molecular structures."
    )

    # Check if the selected drugs have SMILES data
    drugs_with_smiles = [d for d in selected_proteomics_drugs if d in drug_name_to_smiles]
    similar_drugs_with_smiles = [d for d in similar_proteomics_drugs if d in drug_name_to_smiles]

    if drugs_with_smiles and similar_drugs_with_smiles:
        fig = plot_molecular_similarity_heatmap(
            drugs_with_smiles, similar_drugs_with_smiles, drug_name_to_smiles
        )
        st.pyplot(fig)
    else:
        st.warning("No molecular structure data available for the selected drugs.")

    # Display raw similarity data in expander
    with st.expander("View similarity data"):
        # Create a dataframe of all pairwise similarities
        similarity_data = []

        for drug in selected_proteomics_drugs:
            # Get similar drugs for this drug
            drug_data = df_proteomics[df_proteomics["DRUG_NAME_1"] == drug]

            if not drug_data.empty:
                # Sort by similarity score
                drug_data = drug_data.sort_values(
                    by="DCMOA Pearson (r)", ascending=False
                )

                # Add to our collection
                similarity_data.append(drug_data)

        if similarity_data:
            combined_data = pd.concat(similarity_data)
            st.dataframe(combined_data, use_container_width=True)
        else:
            st.write("No similarity data available for the selected drugs.")
