import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D






def plot_data(bar_df, scatter_df, outputImage):
    category_order = ["snp", "snp-complex", "indel", "indel-complex", "large-deletion", "large-insertion", "large-complex"]
    methods_order = ["smooth_10000000_log_clean", "smooth_10000000_log_noClean", "smooth_10000000_noLog_clean", "smooth_10000000_noLog_noClean", "smooth_1_log_clean", "smooth_1_log_noClean", "smooth_1_noLog_clean", "smooth_1_noLog_noClean"]

    bar_df["Category"] = pd.Categorical(bar_df["Category"], categories=category_order, ordered=True)
    bar_df["Method"] = pd.Categorical(bar_df["Method"], categories=methods_order, ordered=True)
    bar_df.sort_values(by=["Category", "Method"], inplace=True)

    scatter_df["Category"] = pd.Categorical(scatter_df["Category"], categories=category_order, ordered=True)
    scatter_df["Method"] = pd.Categorical(scatter_df["Method"], categories=methods_order, ordered=True)
    scatter_df.sort_values(by=["Category", "Method"], inplace=True)

    plt.figure(figsize=(12, 6))

    bar_width = 0.1
    bar_positions = np.arange(len(category_order))

    for idx, method in enumerate(bar_df["Method"].unique()):
        method_data_bar = bar_df[bar_df["Method"] == method]
        plt.bar(bar_positions + idx * bar_width, method_data_bar["fscore"], bar_width, color=f'C{idx}', edgecolor='black', linewidth=0.5)

        method_data_scatter = scatter_df[scatter_df["Method"] == method]
        plt.scatter(bar_positions + idx * bar_width, method_data_scatter["fscore"], color=f'C{idx}', edgecolor='black', marker='o', s=50)

    xticks = []
    for category in category_order:
        numVariants_repeated = int(bar_df[bar_df["Category"] == category]['size'].mean())
        numVariants_nonrepeated = int(scatter_df[scatter_df["Category"] == category]['size'].mean())

        xticks.append(f"{category}\nRepeated: {numVariants_repeated}\nNonRepeated: {numVariants_nonrepeated}")

    plt.xticks(bar_positions + 2 * bar_width, xticks)
    plt.ylabel("F-score")
  #  plt.title(title)

    # Legend for Genotyping tools
    legend_elements_tools = [Patch(facecolor='C' + str(i), edgecolor='black', label=bar_df["Method"].unique()[i]) for i in range(len(bar_df["Method"].unique()))]
    
    # Legend for Region Type
    legend_elements_region = [Patch(facecolor='white', edgecolor='black', label='Bars: Repeated Region'),
                                  Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Circles: Non-Repeated Region')]

    # Create legends
    legend1 = plt.legend(handles=legend_elements_tools, title="Genotyping tools", loc="upper center", bbox_to_anchor=(0.3, -0.13), ncol=len(legend_elements_tools)/2)
    legend2 = plt.legend(handles=legend_elements_region, title="Region Type", loc="upper left", bbox_to_anchor=(0.8, -0.13))
    
    plt.gca().add_artist(legend1)  # We have to add the first legend back manually
    plt.tight_layout()

    # To make sure the legends are not cut off when saving
    plt.savefig(outputImage, dpi=300, bbox_extra_artists=(legend1, legend2), bbox_inches='tight')

if __name__ == "__main__":
    inputFile_repeated = sys.argv[1]
    inputFile_nonrepeated = sys.argv[2]
    outputImage = sys.argv[3]

    bar_df = pd.read_csv(inputFile_repeated)
    scatter_df = pd.read_csv(inputFile_nonrepeated)

    plot_data(bar_df, scatter_df, outputImage)
