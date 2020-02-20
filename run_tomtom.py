"""
This module contains functions to rank motifs using tomtom,
summarize and plot the data.

Requires:
    meme 4.10.0 obtainable from http://meme-suite.org/doc/download.html?man_type=web

Usage:
    python run_tomtom.py <Tf_name> <test_meme_file>  <results_path>
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from past.utils import old_div
import os
import sys
from math import log

import numpy as np
import pandas as pd
from scipy.stats.mstats import winsorize
import seaborn as sns
from matplotlib import pyplot as plt

from .utils import meme_path


def run_tomtom(tf, meme_file, results_folder, figure=False):
    """

    :param tf: Transcription factor name
    :param meme_file: meme input file
    :param results_folder: Link to results folder
    :param figure:
    :return: Saves data to file
    """
    #meme_file = "%s/%s/%s" % (results_folder, tf, tf)
    results_path = "%s/%s" % (results_folder, tf)
    os.system(
        "%s/tomtom -min-overlap 1 -dist %s -evalue -text -thresh 1000 -verbosity 1 %s %s >%s_tomtom.txt" %
        (meme_path, "ed", meme_file, meme_file, results_path))

    clean_tomtom("%s_tomtom.txt" % results_path, "%s.tomtom" % results_path)
    if figure:
        plot_tomtom("%s.tomtom" % results_path, "%s_tomtom.png" % results_path)
        plot_tomtom("%s.tomtom" % results_path, "%s_tomtom.eps" % results_path)


def clean_tomtom(tom_in, tom_out):

    tomtom = pd.read_table(tom_in)

    tomtom["p-value"] += 0.000000000000000000000000000000000000000000000001

    tomtom["p-value"] = -tomtom["p-value"].apply(log)

    tomtom = tomtom[["#Query ID", "Target ID", "p-value"]]

    tomtom.columns = [["Query_ID", "Target_ID", "Score"]]

    tomtom_matrix = tomtom.pivot(index="Query_ID", columns="Target_ID", values="Score")

    #tomtom_matrix.drop_duplicates(inplace=True)

    #tomtom_matrix = tomtom_matrix.T

    #tomtom_matrix.drop_duplicates(inplace=True)

    # Winsorize to reduce effect of outliers
    tomtom_winz = pd.DataFrame(
        winsorize(np.array(tomtom_matrix.values.T.tolist()),
                  limits=0.05), index=tomtom_matrix.index, columns=tomtom_matrix.columns)

    # Normalize
    tomtom_normalized = old_div(tomtom_winz, tomtom_winz.max())

    tomtom_normalized["Average"] = tomtom_normalized.mean()
    tomtom_normalized.sort_values(by="Average", ascending=False, inplace=True)
    tomtom_normalized.loc["Average"] = tomtom_normalized.mean(axis=0)
    tomtom_normalized = tomtom_normalized.T
    tomtom_normalized.sort_values(by="Average", axis=0, ascending=False, inplace=True)
    tomtom_normalized.drop("Average", axis=0, inplace=True)

    # Save the data to file
    tomtom_normalized.to_csv("%s" % tom_out, sep="\t")

    return tomtom_normalized


def plot_tomtom(tom_in, figure_output):

    tomtom_normalized = pd.read_table(tom_in, index_col="Target_ID")

    cg = sns.clustermap(tomtom_normalized, method='single', metric="euclidean",
                        row_cluster=False, linewidths=.15)
    test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    fig = plt.gcf()
    fig.savefig(figure_output, bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    tf = sys.argv[1]
    test_meme_input = sys.argv[2]
    results_path = sys.argv[3]
    figure = True
    run_tomtom(tf, test_meme_input, results_path, figure)
