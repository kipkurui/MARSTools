"""
This module contains functions used to rank motifs using tomtom,
summarize and plot the data.
Requires:
    meme 5.1.1 obtainable from http://meme-suite.org/doc/download.html?man_type=web
Usage:
    python run_tomtom.py <Tf_name> <test_meme_file>  <results_path>
"""

# Import the required modules
import os
import sys
from math import log
import numpy as np
import pandas as pd
from scipy.stats.mstats import winsorize
import seaborn as sns
from matplotlib import pyplot as plt


def run_tomtom(tf, meme_file, results_folder, figure=False):
    """
    :param tf: Transcription factor name
    :param meme_file: meme input file
    :param results_folder: Link to results folder
    :param figure:
    :return: Saves data to file
    """
    # Run Tomtom analysis
    results_path = "%s/%s" % (results_folder, tf)
    os.system(
        "tomtom -min-overlap 1 -dist %s -evalue -text -thresh 1000 -verbosity 1 %s %s > %s_raw.tomtom" %
        ("ed", meme_file, meme_file, results_path))

    # Remove the extra comment lines at the bottom when using meme v 5.1.1
    f = open(results_folder + "/" + tf + "_raw.tomtom", "r")
    tmp = f.readlines()

    # Remove the last few lines which give the meme version and tomtom command used
    with open(results_folder + "/" + tf + "_raw.tomtom", "w") as f:
        for line in tmp:
            if line.startswith("#"):
                continue
            else:
                f.write(line)
    # clean up the results
    clean_tomtom("%s_raw.tomtom" % results_path, "%s.tomtom" % results_path)
    if figure:
        plot_tomtom("%s.tomtom" % results_path, "%s_tomtom.png" % results_path)
        # plot_tomtom("%s.tomtom" % results_path, "%s_tomtom.eps" % results_path)


def clean_tomtom(tom_in, tom_out):
    tomtom = pd.read_csv(tom_in, sep='\t')

    # To avoid zero division, use add the smallest possible value to each

    tomtom["p-value"] += 0.000000000000000000000000000000000000000000000001

    # Get the negative log of the p-value column
    tomtom["p-value"] = -tomtom["p-value"].apply(log)

    # Exatract the important columns to be used for any further analysis
    tomtom = tomtom[["Query_ID", "Target_ID", "p-value"]]

    # Change the column names
    # When renaming columns, use DataFrame.columns = [list], not DataFrame.columns = [[list]]:
    tomtom.columns = ["Query_ID", "Target_ID", "Score"]

    # Pivot the data into a pairwise matrix
    tomtom_matrix = tomtom.pivot(index="Query_ID", columns="Target_ID", values="Score")

    # tomtom_matrix.drop_duplicates(inplace=True)

    # tomtom_matrix = tomtom_matrix.T

    # tomtom_matrix.drop_duplicates(inplace=True)

    # convert the dataframe to an array the winsorize and finally convert back into a dataframe
    # Winsorize to reduce effect of outliers
    tomtom_winz = pd.DataFrame(
        winsorize(np.array(tomtom_matrix.values.T.tolist()),
                  limits=0.05), index=tomtom_matrix.index, columns=tomtom_matrix.columns)

    # Normalize the data
    tomtom_normalized = tomtom_winz.div(tomtom_winz.max())

    tomtom_normalized["Average"] = tomtom_normalized.mean()
    tomtom_normalized.sort_values(by="Average", ascending=False, inplace=True)
    tomtom_normalized.loc["Average"] = tomtom_normalized.mean(axis=0)
    tomtom_normalized = tomtom_normalized.T
    tomtom_normalized.sort_values(by="Average", axis=0, ascending=False, inplace=True)
    tomtom_normalized.drop("Average", axis=0, inplace=True)

    # Save the data to file
    tomtom_normalized.to_csv("%s" % tom_out, sep="\t")

    # return tomtom_normalized


def plot_tomtom(tom_in, figure_output):
    tomtom_normalized = pd.read_csv(tom_in, index_col="Target_ID", sep="\t")

    no_rows, no_cols = tomtom_normalized.shape
    if no_rows <= 30:
        cg = sns.clustermap(tomtom_normalized, method='single', metric="euclidean", row_cluster=False, linewidth=.005,
                            cbar_pos=(0.05, .25, .03, .5), cmap="vlag", standard_scale=1)
    elif no_rows <= 60:
        cg = sns.clustermap(tomtom_normalized, method='single', metric="euclidean", row_cluster=False, linewidth=.005,
                            figsize=(13, 13), cbar_pos=(0.05, .25, .03, .5), cmap="vlag", standard_scale=1)
    else:
        cg = sns.clustermap(tomtom_normalized, method='single', metric="euclidean", row_cluster=False, linewidth=.005,
                            figsize=(16, 16), cbar_pos=(0.05, .25, .03, .5), cmap="vlag", standard_scale=1)
    test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    fig = plt.gcf()
    fig.savefig(figure_output, bbox_inches='tight', dpi=100)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    tf = sys.argv[1]
    test_meme_input = sys.argv[2]
    results_path = sys.argv[3]
    figure = True
    run_tomtom(tf, test_meme_input, results_path, figure)
