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

D_NAME = "DRUG_NAME"
G_DRIVER = "Driver_Gene"
G_NAME = "GENE"
D_ACT = "ACT_VALUE"


def radar_plot(df):
    """
    Create radar plots from drug activity data.

    Parameters:
    -----------
    df_gene_drug_filtered : pandas DataFrame
        DataFrame with columns DRUG_NAME, GENE, ACT_VALUE, ACT_TYPE
    output_file : str, optional
        Path to save the plot. If None, the plot will be displayed.
    """
    # Group by ACT_TYPE
    act_types = df["ACT_TYPE"].unique()

    # Set up the figure
    fig = plt.figure(figsize=(15, 10))

    # Create a radar chart for each ACT_TYPE
    for i, act_type in enumerate(act_types):
        # Filter data for this activity type
        act_data = df[df["ACT_TYPE"] == act_type]

        # If there are duplicate DRUG_NAME entries, take the mean of ACT_VALUE
        act_data = act_data.groupby("DRUG_NAME")["ACT_VALUE"].mean().reset_index()

        # Get the drug names and activity values
        drug_names = act_data["DRUG_NAME"].tolist()
        values = act_data["ACT_VALUE"].tolist()

        # Create subplot
        ax = fig.add_subplot(1, len(act_types), i + 1, polar=True)

        # Number of variables
        N = len(drug_names)

        # If we don't have any data, skip this subplot
        if N == 0:
            continue

        # What will be the angle of each axis in the plot
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop

        # Draw one axis per variable and add labels
        plt.xticks(angles[:-1], drug_names, color="black", size=8)

        # Draw the activity values
        values += values[:1]  # Close the loop
        ax.plot(angles, values, linewidth=2, linestyle="solid", label=act_type)
        ax.fill(angles, values, alpha=0.25)

        # Add legend
        ax.legend(loc="upper right", bbox_to_anchor=(0.1, 0.1))

        # Set chart title
        plt.title(f"Activity Type: {act_type}", size=12, color="blue", y=1.1)

    # Adjust layout
    plt.tight_layout()

    return fig


def plot_graph(df_drug_gene):
    # Import necessary libraries

    # Create a bipartite graph
    B = nx.Graph()

    # Extract unique drugs and genes
    drugs = df_drug_gene[D_NAME].unique()
    genes = df_drug_gene[G_NAME].unique()

    # Add nodes with the bipartite attribute
    B.add_nodes_from(drugs, bipartite=0)  # Drugs are group 0
    B.add_nodes_from(genes, bipartite=1)  # Genes are group 1

    # Add edges from drug-gene pairs in the dataframe
    edges = list(zip(df_drug_gene[D_NAME], df_drug_gene[G_NAME]))
    B.add_edges_from(edges)

    # Check the size of the graph
    print(f"Graph has {B.number_of_nodes()} nodes and {B.number_of_edges()} edges")

    # If the graph is too large, consider filtering
    if B.number_of_nodes() > 100:
        # Get the top N drugs by degree (number of targets)
        drug_degrees = [(node, d) for node, d in B.degree() if node in drugs]
        drug_degrees.sort(key=lambda x: x[1], reverse=True)
        top_drugs = [node for node, d in drug_degrees[:15]]  # Get top 15 drugs

        # Create a subgraph with only the top drugs and their targets
        sub_nodes = set(top_drugs)
        for drug in top_drugs:
            sub_nodes.update(B.neighbors(drug))

        B = B.subgraph(sub_nodes)
        print(
            f"Filtered graph has {B.number_of_nodes()} nodes and {B.number_of_edges()} edges"
        )

    # Create the layout
    pos = nx.drawing.layout.bipartite_layout(B, drugs)

    # Plot the graph with custom settings
    fig = plt.figure(figsize=(16, 10))

    # Draw drugs (group 0)
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=list(set(drugs) & set(B.nodes)),
        node_color="lightblue",
        node_size=200,
        alpha=0.8,
    )

    # Draw genes (group 1)
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=list(set(genes) & set(B.nodes)),
        node_color="lightgreen",
        node_size=100,
        alpha=0.8,
    )

    # Draw edges
    nx.draw_networkx_edges(B, pos, width=0.5, alpha=0.5, edge_color="gray")

    # Add labels with smaller font for readability
    nx.draw_networkx_labels(B, pos, font_size=8)

    plt.title("Drug-Gene Network")
    plt.tight_layout()
    plt.axis("off")
    # plt.show()
    return fig


PROCESSED_DATA_PATH = Path("data/DrugCentral/processed_data/")

df_cancer_gene = pd.read_csv(PROCESSED_DATA_PATH / "gene.id.csv")
df_gene_drug = pd.read_table(
    PROCESSED_DATA_PATH / "drug.target.interaction.fda.cosmic.cancer.type.tsv"
)
df_gene_drug = df_gene_drug.query("ACT_TYPE in ['Kd', 'IC50']")

drug_name_id_map = pd.read_csv(PROCESSED_DATA_PATH / "drug.structure.id.csv")
cancers_options = df_cancer_gene["CANCER"].unique()
cancer = st.selectbox("Select cancer", cancers_options)

genes_in_cancer = df_cancer_gene.query(f"CANCER=='{cancer}'")["GENE"]

st.write(f"Genes number: {genes_in_cancer.nunique()}")
# st.write(df)
# genes
df_gene_drug_filtered = df_gene_drug.query("GENE in @genes_in_cancer")
df_gene_drug_filtered

# df_gene_drug_filtered have columns DRUG_NAME GENE ACT_VALUE ACT_TYPE
#
fig = radar_plot(df_gene_drug_filtered)
fig

fig = plot_graph(df_gene_drug_filtered)
fig
