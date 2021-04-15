import os
import sys
import shutil

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from .utils import tab2fasta, rotate_image, meme_path


def run_gimme(tf, user_motif, chip_seq_list, results_path, genome, figure=False):
    """
    """
    files_path = '%s/%s' % (results_path, tf)
    gimme_in = "%s/%s_gimme_metrics.txt" % (results_path, tf)
    gimme_out = "%s/%s.gimme" % (results_path, tf)

    # Remove the metrics file if it exists
    if os.path.exists("%s/%s_gimme_metrics.txt" % (results_path, tf)):
        os.remove("%s/%s_gimme_metrics.txt" % (results_path, tf))
    # Loop through the chip peaks getting the metrics    
    for chip_seq in chip_seq_list:
        file_name = chip_seq.split('/')[-1].split('.')[0]

        # Extract the positive and negative sequence
        tab2fasta(chip_seq, '%s/%s.fa' % (results_path, file_name), '%s/%s.bg' % (results_path, file_name))

        # Convert the user motif into a pfm format
        gimme_mot = '%s/%s.pfm' % (results_path, tf)
        meme2gimme(user_motif, gimme_mot)

        # Run the gimme enrichment analysis
        os.system("gimme motifs %s/%s.fa %s/%s -b %s/%s.bg -p %s/%s.pfm --known -g %s" %
                  (results_path, file_name, results_path, tf, results_path, file_name, results_path, tf, genome))
        # Append all the chip results into one file
        os.system("cat %s/%s/gimme.roc.report.txt >> %s/%s_gimme_metrics.txt" % (results_path, tf, results_path, tf))

        # Clean up the files
        os.remove("%s/%s.fa" % (results_path, file_name))
        os.remove("%s/%s.bg" % (results_path, file_name))
        shutil.rmtree("%s/%s" % (results_path, tf))
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
            outs.write(
                "Motif\t# matches\t# matches background\tP-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% "
                "FPR\tRecall at 10% FDR\n")
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
    new_gimme = new_gimme.sort_values(by="ROC AUC", axis=0, ascending=False)
    labels = new_gimme.index
    x = 10
    if len(labels) > 50:
        x = 15
    elif len(labels) < 10:
        x = 5

    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(x, 10), sharex='none', dpi=100)
    a = sns.barplot(x=labels, y=new_gimme["ROC AUC"],
                    palette='colorblind', order=labels, ax=ax1)
    b = sns.barplot(x=labels, y=new_gimme["PR AUC"],
                    palette="colorblind", order=labels, ax=ax2)
    c = sns.barplot(x=labels, y=new_gimme["Enr. at 1% FPR"],
                    palette="colorblind", order=labels, ax=ax3)
    d = sns.barplot(x=labels, y=new_gimme["Recall at 10% FDR"],
                    palette="colorblind", order=labels, ax=ax4)
    e = sns.barplot(x=labels, y=new_gimme["P-value"],
                    palette="colorblind", order=labels, ax=ax5)

    a.set_xlabel("")
    b.set_xlabel("")
    c.set_xlabel("")
    d.set_xlabel("")
    e.set_xticklabels(labels, rotation=90)
    e.set_xlabel("Motif", fontdict={'fontsize': 12, 'fontweight': 'semibold'})
    sns.despine()
    f.savefig(figure_out, bbox_inches='tight')


def meme2gimme(meme, gimme):
    with open(meme) as motif:
        with open(gimme, 'w') as gmotif:
            for line in motif:
                if line.startswith("MOTIF"):
                    if len(line.split(" ")) > 2:
                        gmotif.write(">" + line.split(" ")[1] + "\n")
                    else:
                        gmotif.write(">" + line.split(" ")[1])
                elif line.startswith('letter-probability'):
                    continue
                elif line.startswith('  '):
                    a = line.split()
                    if len(a) > 0:
                        gmotif.write(a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\n')
                    else:
                        continue
                elif line.startswith("\n"):
                    continue
                else:
                    continue


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    TF = sys.argv[1]  # type: str
    CHIP_SEQ_LIST = sys.argv[2]
    MEME_INPUT = sys.argv[3]
    RESULTS_PATH = sys.argv[4]
    GENOME = sys.argv[5]
    run_gimme(TF, MEME_INPUT, CHIP_SEQ_LIST, RESULTS_PATH, GENOME)