import os
import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from utils import tab2fasta, rotate_image, meme_path


def run_gimme(tf, user_motif, chip_seq_list, results_path, figure=False):
    files_path = '%s/%s' % (results_path, tf)
    gimme_in = "%s/%s_gimme_metrics.txt" % (results_path, tf)
    gimme_out = "%s/%s.gimme" % (results_path, tf)
    for chip_seq in chip_seq_list:
        file_name = chip_seq.split('/')[-1].split('.')[0]
        tab2fasta(chip_seq, '%s/%s.fa' % (results_path, file_name), '%s/%s.bg' % (results_path, file_name))
        gimme_mot = '%s/%s.motif' % (results_path, tf)
        meme2gimme(user_motif, gimme_mot)

        os.system("%s/gimme roc %s %s/%s.fa %s/%s.bg >>%s/%s_gimme_metrics.txt" %
                  (meme_path, gimme_mot, results_path, file_name, results_path, file_name, results_path, tf))
    # import glob
    # for i in glob.glob('tmp/*'):
    #     os.remove(i)

    if figure:
        plot_histogram_gimme(gimme_in, gimme_out, "%s/%s_gimme.png" % (results_path, tf))
        plot_histogram_gimme(gimme_in, gimme_out, "%s/%s_gimme.eps" % (results_path, tf))

        rotate_image("%s/%s_gimme.png" % (results_path, tf), "%s/%s_gimme_rot.png" % (results_path, tf))


def clean_gimme(gimme_in, gimme_out):
    """
    Function that cleans the GIMME output for the purpose of plotting the data
    """
    with open(gimme_in) as gim:
        with open(gimme_out, 'w') as outs:
            outs.write("Motif	ROC AUC	MNCP	Enr. at 5% FDR	Max enr.	Recall at 10% FDR\n")
            for line in gim:
                if line.startswith("Motif"):
                    continue
                else:
                    outs.write(line)


def plot_histogram_gimme(gimme_in, gimme_out, figure_out):
    """
    This function plots all the assess motifs output on a single plot and
    sorts the data based on a statistical function used.

    The function could be modified for flexibility in that it can be used
    to generate both plotly and seaborn figures. Also, a choice of the
    statistical tool used for sorting the data can be added through the
    parameters.
    """

    clean_gimme(gimme_in, gimme_out)
    gimme = pd.read_table(gimme_out)
    new_gimme = pd.pivot_table(gimme, index=['Motif'])
    new_gimme = new_gimme.sort(columns="ROC AUC", axis=0, ascending=False)
    labels = new_gimme.index
    x = 10
    if len(labels) > 50:
        x = 15
    elif len(labels) < 10:
        x = 5
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(x, 10), sharex=True)
    a = sns.barplot(x=labels, y=new_gimme["ROC AUC"],
                    palette='colorblind', x_order=labels, ax=ax1)
    b = sns.barplot(x=labels, y=new_gimme["MNCP"],
                    palette="colorblind", x_order=labels, ax=ax2)
    c = sns.barplot(x=labels, y=new_gimme["Enr. at 5% FDR"],
                    palette="colorblind", x_order=labels, ax=ax3)
    d = sns.barplot(x=labels, y=new_gimme["Max enr."],
                    palette="colorblind", x_order=labels, ax=ax4)
    e = sns.barplot(x=labels, y=new_gimme["Recall at 10% FDR"],
                    palette="colorblind", x_order=labels, ax=ax4)
    d.set_xticklabels(labels, rotation=90)
    sns.despine()
    f.savefig(figure_out, bbox_inches='tight')


def meme2gimme(meme, gimme):
    with open(meme) as motif:
        with open(gimme, 'w') as gmotif:
            for line in motif:
                if line.startswith("MOTIF"):
                    if len(line.split(" ")) > 2:
                        gmotif.write(">"+line.split(" ")[1]+"\n")
                    else:
                        gmotif.write(">"+line.split(" ")[1])
                elif line.startswith('letter-probability'):
                    continue
                elif line.startswith('  '):
                    a = line.split()
                    if len(a) > 0:
                        gmotif.write(a[0]+'\t'+a[1]+'\t'+a[2]+'\t'+a[3]+'\n')
                    else:
                        continue
                elif line.startswith("\n"):
                    continue
                else:
                    continue

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print __doc__
        sys.exit(1)
    tf = sys.argv[1]
    chip_seq_list = sys.argv[2]
    test_meme_input = sys.argv[3]
    results_path = sys.argv[4]

    run_gimme(tf, test_meme_input, chip_seq_list, results_path)
