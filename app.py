from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np


from plots import (
    plot_graph,
    plot_molecular_similarity_heatmap,
    plot_proteomics_barplots,
    radar_plot,
)

D_NAME = "DRUG_NAME"
G_DRIVER = "Driver_Gene"
G_NAME = "GENE"
D_ACT = "ACT_VALUE"


PROCESSED_DATA_PATH = Path("data/DrugCentral/processed_data/")

df_cancer_gene = pd.read_csv(PROCESSED_DATA_PATH / "gene.id.csv")
df_gene_drug = pd.read_table(
    PROCESSED_DATA_PATH / "drug.target.interaction.fda.cosmic.cancer.type.tsv"
)
df_proteomics = pd.read_csv(
    "data/DeepCoverageMOA/analyzed_data/deepcoveragemoa.lincs.csv"
)

df_drug_to_smiles = pd.read_table("data/DrugCentral/raw_data/structures.smiles.tsv")[
    ["INN", "SMILES"]
]

drug_name_to_smiles = {t.INN: t.SMILES for t in df_drug_to_smiles.itertuples()}


df_gene_drug = df_gene_drug.query("ACT_TYPE in ['Kd', 'IC50']")

drug_name_id_map = pd.read_csv(PROCESSED_DATA_PATH / "drug.structure.id.csv")
cancers_options = df_cancer_gene["CANCER"].unique()
cancer = st.selectbox("Select cancer", cancers_options)

genes_in_cancer = df_cancer_gene.query(f"CANCER=='{cancer}'")["GENE"]

st.write(f"Genes number: {genes_in_cancer.nunique()}")
# st.write(df)
# genes
df_gene_drug_filtered = df_gene_drug.query("GENE in @genes_in_cancer")

# allow user to select row from df_gene_drug_filtered and save selected drug names to variable


# DRUG-GENE plots
fig = radar_plot(df_gene_drug_filtered)
fig

fig = plot_graph(df_gene_drug_filtered)
fig


st.write("### Select drugs of interest")

# Create a multi-select widget to select drugs
unique_drugs = df_gene_drug_filtered[D_NAME].unique()
drugs_with_proteomics = set(df_proteomics["DRUG_NAME_1"].unique()) & set(
    df_gene_drug_filtered[D_NAME].unique()
)

# unique_drugs = [drug for drug in unique_drugs if drug in df_proteomics["DRUG_NAME_1"]]


selected_drugs = st.multiselect("Select drugs", drugs_with_proteomics)


if selected_drugs:
    figs, similar_drugs = plot_proteomics_barplots(df_proteomics, selected_drugs)
    for fig in figs:
        fig
    # st.write(df_proteomics)

    fig = plot_molecular_similarity_heatmap(
        selected_drugs, similar_drugs, drug_name_to_smiles
    )
    fig