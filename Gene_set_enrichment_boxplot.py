import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import pandas as pd

def analyze_and_plot(data, cluster_col, score_col, list_colors, save_path):
    """
    This function performs a Mann-Whitney U test between groups in the data and plots the results along with box plots.

    Parameters:
    data (DataFrame): The dataset containing the scores and cluster information.
    cluster_col (str): The name of the column in 'data' that contains cluster identifiers.
    score_col (str): The name of the column in 'data' that contains the scores to analyze.
    list_colors (list): The list of colors to use for the box plots.
    save_path (str): The path where the plot image will be saved.

    Returns:
    None: The function saves the plot image to 'save_path' and displays the plot.
    """

    # Prepare a dictionary to store the results of the tests
    test_results = {}

    # Identify the unique groups and count of samples per group
    unique_groups = data[cluster_col].unique()

    # Data for the comparison group (group 3)
    group3_data = data[data[cluster_col] == "3"][score_col]

    # Perform the Mann-Whitney U test between group 3 and all other groups
    for group in unique_groups:
        if group == "3":
            continue  # Skip the comparison with itself
        
        # Data for the current group
        current_group_data = data[data[cluster_col] == group][score_col]
        
        # Perform the test
        stat, p_value = mannwhitneyu(group3_data, current_group_data, alternative='two-sided')
        
        # Save the results
        test_results[group] = {'U_stat': stat, 'p_value': p_value}

    # Create the box plot with horizontal comparison lines
    plt.figure(figsize=(12, 7))

    # Box plot
    sns.boxplot(x=cluster_col, y=score_col, data=data, palette=list_colors)  # Change the color palette if necessary

    # Add lines and annotations for statistical comparisons
    y_min, y_max = plt.ylim()  # Get the current y limits
    y_diff = y_max - y_min  # Calculate the difference

    for i, (group, stats) in enumerate(test_results.items()):
        group = int(group)
        # y coordinate for the line, increases with each iteration to prevent line overlap
        y_line = y_max + (i + 1) * (y_diff * 0.1)
        
        # Draw the horizontal line from the current group to group 3
        plt.hlines(y_line, group, 3, color='black', linestyle='--')
        
        # Position for the p-value text
        x_text = (group + 3) / 2
        plt.text(x_text, y_line, f"p={stats['p_value']:.2e}", ha='center', va='bottom')

    # Update the upper y limit to make room for lines and annotations
    plt.ylim(top=y_line + (y_diff * 0.1))

    # Labels and title
    plt.title('Distribution of gene set scores per cluster with statistical comparisons')
    plt.ylabel('Gene set scores')
    plt.xlabel('Leiden Cluster')

    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()
