import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

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


def plot_graph(df_drug_gene, node_spacing=0.3):
    """
    Create a bipartite graph visualization of drug-gene interactions with improved spacing.

    Parameters:
    -----------
    df_drug_gene : pandas DataFrame
        DataFrame containing drug-gene interactions
    node_spacing : float, default 0.3
        Vertical spacing between nodes in the layout

    Returns:
    --------
    matplotlib figure
    """
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
    st.write(f"Graph has {B.number_of_nodes()} nodes and {B.number_of_edges()} edges")

    # Make sure genes contains only actual genes (not drugs)
    genes = [node for node in genes if node in B.nodes()]

    # Create a custom layout that spaces nodes better
    pos = {}

    # Position drugs on the left side with improved vertical spacing
    drug_count = len(drugs)
    for i, drug in enumerate(drugs):
        # Distribute drugs evenly along the y-axis
        y_pos = (i - drug_count / 2) * node_spacing
        pos[drug] = np.array([-1.0, y_pos])

    # Position genes on the right side
    gene_count = len([g for g in genes if g in B.nodes()])
    gene_list = [g for g in genes if g in B.nodes()]

    for i, gene in enumerate(gene_list):
        # Distribute genes evenly along the y-axis
        y_pos = (i - gene_count / 2) * node_spacing
        pos[gene] = np.array([1.0, y_pos])

    # Plot the graph with custom settings
    fig = plt.figure(figsize=(16, max(10, drug_count * 0.6)))

    # Draw drugs (group 0)
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=[d for d in drugs if d in B.nodes()],
        node_color="lightblue",
        node_size=300,
        alpha=0.8,
    )

    # Draw genes (group 1)
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=[g for g in genes if g in B.nodes()],
        node_color="lightgreen",
        node_size=150,
        alpha=0.8,
    )

    # Draw edges with curved lines to reduce overlap
    nx.draw_networkx_edges(
        B,
        pos,
        width=0.8,
        alpha=0.5,
        edge_color="gray",
        connectionstyle="arc3,rad=0.1",  # Add a slight curve to edges
    )

    # Add labels with adjusted font sizes and positions
    # Drug labels
    drug_labels = {node: node for node in drugs if node in B.nodes()}
    nx.draw_networkx_labels(
        B,
        {
            k: (v[0] - 0.1, v[1]) for k, v in pos.items() if k in drugs
        },  # Slightly offset labels
        labels=drug_labels,
        font_size=10,
        font_weight="bold",
        horizontalalignment="right",
    )

    # Gene labels
    gene_labels = {node: node for node in genes if node in B.nodes()}
    nx.draw_networkx_labels(
        B,
        {
            k: (v[0] + 0.1, v[1]) for k, v in pos.items() if k in genes
        },  # Slightly offset labels
        labels=gene_labels,
        font_size=8,
        horizontalalignment="left",
    )

    plt.title("Drug-Gene Network", fontsize=16)
    plt.tight_layout()
    plt.axis("off")

    return fig


def plot_proteomics_barplots(df_proteomics, selected_drugs, top_n=10):
    """
    Create barplots for each selected drug showing the top similar drugs from proteomics data.

    Parameters:
    -----------
    df_proteomics : pandas DataFrame
        DataFrame containing proteomics data with similarity scores
    selected_drugs : list
        List of drug names to plot
    top_n : int, default 10
        Number of top similar drugs to display

    Returns:
    --------
    list of matplotlib figures
    """
    figures = []
    similar_drugs = []
    # For each selected drug
    for drug in selected_drugs:
        # Filter the proteomics data for the selected drug
        drug_data = df_proteomics[df_proteomics["DRUG_NAME_1"] == drug]

        # If no data found for this drug, skip
        if drug_data.empty:
            print(f"No proteomics data found for {drug}")
            continue

        # Sort by similarity score (assuming it's in a column named 'SIMILARITY_SCORE')
        # Adjust column name if needed
        score_column = "DCMOA Pearson (r)"

        # Sort and get top N similar drugs
        drug_data = drug_data.sort_values(by=score_column, ascending=False).head(top_n)

        # Create a barplot
        fig, ax = plt.subplots(figsize=(12, 6))

        # Determine which column contains the similar drug names
        drug_name_col = "DRUG_NAME_2"

        # Create the barplot
        bars = ax.barh(drug_data[drug_name_col], drug_data[score_column])

        # Add labels and title
        ax.set_xlabel("Similarity Score")
        ax.set_ylabel("Similar Drugs")
        ax.set_title(f"Top {top_n} Drugs Similar to {drug}")

        # Add values on bars
        for bar in bars:
            width = bar.get_width()
            ax.text(
                width + 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{width:.3f}",
                ha="left",
                va="center",
            )

        # Adjust layout
        plt.tight_layout()
        figures.append(fig)
        similar_drugs.extend(drug_data[drug_name_col])

    return figures, list(set(similar_drugs))


def calculate_similarity_score(mol1, mol2):
    """Calculate Tanimoto similarity between two molecules."""
    return DataStructs.TanimotoSimilarity(mol1, mol2)


def plot_molecular_similarity_heatmap(selected_drugs, similar_drugs, name_to_smiles):
    """
    Create a heatmap showing molecular similarity between drugs.

    Parameters:
    -----------
    selected_drugs : list
        List of selected drug names
    similar_drugs : list
        List of similar drug names from proteomics analysis
    name_to_smiles : dict
        Dictionary mapping drug names to SMILES strings

    Returns:
    --------
    matplotlib figure
    """

    def convert_smiles_to_molecules_and_calculate_score(smiles_pair):
        assert len(smiles_pair) == 2

        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_pair]

        fpgen = AllChem.GetRDKitFPGenerator()
        fps = [fpgen.GetFingerprint(x) for x in mol_list]

        similarity = calculate_similarity_score(fps[0], fps[1])

        return similarity

    # Create a dictionary of drug name to RDKit molecule
    combined_drugs = list(set(selected_drugs) | set(similar_drugs))

    # Calculate similarity matrix
    n = len(combined_drugs)
    similarity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                similarity_matrix[i, j] = 1.0  # Self-similarity is 1
            else:
                if (
                    combined_drugs[i] in name_to_smiles
                    and combined_drugs[j] in name_to_smiles
                ):
                    fp1 = name_to_smiles[combined_drugs[i]]
                    fp2 = name_to_smiles[combined_drugs[j]]
                    similarity_matrix[i, j] = (
                        convert_smiles_to_molecules_and_calculate_score([fp1, fp2])
                    )

    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(similarity_matrix, cmap="viridis")

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Tanimoto Similarity", rotation=-90, va="bottom")

    # Add drug names as ticks
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(combined_drugs, rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(combined_drugs, fontsize=10)

    # Add similarity values in the cells
    for i in range(n):
        for j in range(n):
            text_color = "white" if similarity_matrix[i, j] < 0.7 else "black"
            ax.text(
                j,
                i,
                f"{similarity_matrix[i, j]:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=9,
            )

    ax.set_title("Structural Similarity (Tanimoto) Between Drugs", fontsize=14)
    plt.tight_layout()

    return fig
