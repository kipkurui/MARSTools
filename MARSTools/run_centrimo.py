"""
Author: Caleb Kibet

run_centrimo.py contains functions to run a motif enrichment analysis using CentriMo,
summarize and plot the data.

Requires:
    meme 4.10.0 obtainable from http://meme-suite.org/doc/download.html?man_type=web

Takes as input:
    TF name
    A list of ChIP-seq files formatted in TAB rather than FASTA format.
     Provide a folder with the files
    A motif file in MEME format
    A repository to put results

Usage:
    python run_centrimo.py <Tf_name> <chip-seq_list> <test_meme_file> >results_path>
    eg: python run_centrimo.py Cjun  <test_meme_file> >results_path>
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import


from builtins import map
from builtins import range

import os
import sys
import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .utils import tab2fasta, mkdir_p, meme_path, BASE_DIR


def run_centrimo(tf_name, chip_seq_list, test_meme_input, files_path, figure=False):
    """
    Given a list paths to test data and a MEME file, perform motif enrichment analysis
    of the motifs in the sequences and use the output to rank the motif performance

    Assumes that half the input data are test and the rest are used as background

    :param tf_name:
    :param chip_seq_list:
    :param test_meme_input:
    :param files_path:
    :param figure:
    :return:
    """
    # If chip_seq list is greater than 10, randomly sample 10.
    import random
    if len(chip_seq_list) > 10:
        random.seed(10)
        chip_seq_list = random.sample(chip_seq_list, 10)


    test_list = [["Motif"]]  # get list of motifs in meme file
    with open(test_meme_input) as motifs:
        for motif_line in motifs:
            if motif_line.startswith("MOTIF"):
                test_list.append([motif_line.split()[1]])

    def get_centrimo_list(chip_name):
        """
        Extracts important details from a CentriMo run output
        """
        with open(chip_name) as cent:
            temp_dict = {}
            for line in cent:
                if line.startswith("#"):
                    continue
                else:
                    temp_dict[line.split()[1]] = float(line.split()[5])*-1

        for mot in range(len(test_list)):
            if test_list[mot][0] in temp_dict:
                test_list[mot].append(temp_dict[test_list[mot][0]])
            elif test_list[mot][0] == "Motif":
                continue
            else:
                test_list[mot].append(0)

    for chip_seq in chip_seq_list:
        # Get file name from the path
        file_name = chip_seq.split('/')[-1].split('.')[0]

        test_list[0].append(file_name)
        tmp_path = '%s/%s/tmp' % (files_path, tf_name)
        mkdir_p(tmp_path)

        # Convert the tab file to FASTA
        # TODO: Also allow input as FASTA
        tab2fasta(chip_seq, '%s/%s.fa' % (tmp_path, file_name), '%s/%s.bg' % (tmp_path, file_name))

        os.system("%s/fasta-get-markov  %s/%s.fa %s/%s.fa.bg" % (meme_path, tmp_path, file_name, tmp_path, file_name))
        os.system("%s/centrimo --oc %s/%s --verbosity 1 --local --optimize_score --score 5.0 "
                  "--ethresh 100000.0 --neg %s/%s.bg --bfile %s/%s.fa.bg %s/%s.fa %s" %
                  (meme_path, tmp_path, file_name, tmp_path, file_name, tmp_path,
                   file_name, tmp_path, file_name, test_meme_input))

        get_centrimo_list("%s/%s/centrimo.txt" % (tmp_path, file_name))

        # Delete all the temporary files
        import shutil
        shutil.rmtree('%s/%s/' % (tmp_path, file_name))

        for i in glob.glob('%s/*' % tmp_path):
            os.remove(i)
        os.rmdir(tmp_path)

    test_list[0].append("Average")
    for i in range(1, len(test_list)):
        test_list[i].append(np.mean(test_list[i][1:]))
    test_list.sort(key=lambda x: x[-1], reverse=True)

    with open('%s/%s_centrimo.txt' % (files_path, tf_name), 'w') as cent_out:
        for i in test_list:
            cent_out.writelines('\t'.join(map(str, i)) + '\n')
    centrimo_df = pd.read_table('%s/%s_centrimo.txt' % (files_path, tf_name), index_col=0)
    del centrimo_df['Average']

    centrimo_normalized = centrimo_df / centrimo_df.max()
    centrimo_normalized = centrimo_normalized.replace(to_replace='NaN', value=0)
    average_column = centrimo_normalized.T.mean()
    average_column = average_column.to_frame(name="Average")

    # Add the average column to DataFramne
    centrimo_normalized = centrimo_normalized.T.append(average_column.T).T

    centrimo_normalized.sort_values(by="Average", axis=0, ascending=False, inplace=True)
    centrimo_normalized.to_csv('%s/%s_centrimo_norm.txt' % (files_path, tf_name), sep="\t")

    cent_path = '%s/%s_centrimo_norm.txt' % (files_path, tf_name)
    if figure:
        plot_centrimo(cent_path, '%s/%s_centrimo.png' % (files_path, tf_name))
        plot_centrimo(cent_path, '%s/%s_centrimo.eps' % (files_path, tf_name))


def plot_centrimo(centrimo_in, figure_output):
    """

    :param centrimo_in:
    :param figure_output:
    :return:
    """
    centrimo_table = pd.read_table(centrimo_in, index_col=0)
    centrimo_table.sort_values(by="Average", axis=0, ascending=False, inplace=True)

    cg = sns.clustermap(centrimo_table, method='single', metric="euclidean",
                        row_cluster=False, linewidths=.15)
    test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    f = plt.gcf()
    f.savefig(figure_output, bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    TF = sys.argv[1]
    CHIP_SEQ_LIST = sys.argv[2]
    TEST_MEME = sys.argv[3]
    RESULTS_PATH = sys.argv[4]
    run_centrimo(TF, CHIP_SEQ_LIST, TEST_MEME, RESULTS_PATH)