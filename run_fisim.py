"""
Author: Caleb Kibet

run_fisim.py contains functions to evaluate and rank motifs based on consistency using
FISim tools, summarize and plot the results.

Requires:
    meme 4.10.0 obtainable from http://meme-suite.org/doc/download.html?man_type=web

Takes as input:
    TF name
    A motif file in MEME format
    A repository to output results

Usage:
    python run_centrimo.py <Tf_name> <chip-seq_list> <test_meme_file> >results_path>
    eg: python run_centrimo.py Cjun  <test_meme_file> >results_path>
"""

import os
import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from utils import BASE_DIR


def run_fisim(tf, meme_file, results_folder, figure=False):
    """
    Run FISim motif comparison and
    :param meme_file:
    :param tf:
    :param results_folder:
    :param figure:
    :return:
    """

    tf_path = "%s/%s" % (results_folder, tf)
    os.chdir("%s/MARSTools/FISIM" % BASE_DIR)  # May need to rethink the location of the scripts
    os.system("python fisim.py -fileList %s -o %s.fisim -core -ID" % (meme_file, tf_path))
    os.chdir(BASE_DIR)

    fisim_df = pd.read_table(tf_path+'.fisim', index_col="Motifs")
    # Eliminate Duplicates
    fisim_df = fisim_df.drop_duplicates()
    fisim_df = fisim_df.T.drop_duplicates().T
    fisim_df.sort_values("Average", axis=0, ascending=False, inplace=True)

    fisim_df.to_csv(tf_path + ".fisim", sep="\t")
    if figure:
        plot_heatmap_fisim(tf_path + '.fisim', tf_path + '_fisim.png')
        plot_heatmap_fisim(tf_path + '.fisim', tf_path + '_fisim.eps')


def plot_heatmap_fisim(fisim_input, figure_output):
    """
    Generate a figure that depict teh ranking and the clustering of the FISIM motifs
    :param fisim:
    :return:
    """
    fisim_df = pd.read_table(fisim_input, index_col="Motifs")
    # Eliminate Duplicates
    fisim_df = fisim_df.drop_duplicates()
    fisim_df = fisim_df.T.drop_duplicates().T
    fisim_df.sort_values("Average", axis=0, ascending=False, inplace=True)
    fisim_df = fisim_df.T.ix[:-1].T
    cg = sns.clustermap(fisim_df, method='single', metric="euclidean", row_cluster=False, linewidths=.15)
    test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    fig = plt.gcf()
    fig.savefig(figure_output, bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print __doc__
        sys.exit(1)
    tf = sys.argv[1]
    test_meme_input = sys.argv[2]
    results_path = sys.argv[3]
    figure = True
    run_fisim(tf, test_meme_input, results_path, figure)
