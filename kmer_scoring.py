"""
Author: Caleb Kibet

A module for scoring sequences using k-mers used by Assess_by_score

"""

from builtins import range
import pandas as pd

from MARSTools.utils import revcompl


def get_kmer_dict_rev(kmerscore, kmer_name):
    """
    Given a PBM enrichment scores file, create a k-mer dictionary
    of both forward and reverse sequences

    :param kmerscore:
    :param kmer_name:
    :return:
    """
    kmer_df_fw = pd.read_table(kmerscore, index_col="8-mer", usecols=["8-mer", "E-score"])
    kmer_df_fw.fillna(0, inplace=True)
    kmer_df_rv = pd.read_table(kmerscore, index_col="8-mer.1", usecols=["8-mer.1", "E-score"])
    kmer_df_rv.index.name = "8-mer"
    kmer_df_rv.fillna(0, inplace=True)
    combined_kmers = kmer_df_fw.append(kmer_df_rv)
    combined_kmers_dict = combined_kmers.to_dict()["E-score"]

    return combined_kmers_dict, kmer_name


def score_kmer(kmer, kmerdict):
    """
    Simple function to score sequences given a k-mer dictionary
    of forward k-mers

    :param kmer:
    :param kmerdict:
    :return:
    """
    score = 0
    if kmer in kmerdict:
        score = float(kmerdict[kmer])
    else:
        kmer2 = revcompl(kmer)
        score = float(kmerdict[kmer2])
    return score


def sum_kmer_score(kmerdict, seq):
    """
    Given a k-mer dictionary and a sequence, calculate sum occupancy

    :param kmerdict: A dictionary of k-mers
    :param seq: Sequence
    :return: Sequence occupancy score
    """
    k_mers = find_kmers(seq, 8)
    total_seq_score = 0
    for kmer in k_mers:
        if kmer in kmerdict:
            kmer_score = float(kmerdict[kmer])
        else:
            kmer2 = revcompl(kmer)
            kmer_score = float(kmerdict[kmer2])
        total_seq_score += kmer_score
    return total_seq_score


def max_kmer_score(kmerdict, seq):
    """

    Given a sequence and a dictionary of k-mer scores, calculate
    maximum sequence occupancy

    :param kmerdict:
    :param seq:
    :return:
    """
    k_mers = find_kmers(seq, 8)
    kmer_scores_list = []
    for kmer in k_mers:
        if kmer in kmerdict:
            score = float(kmerdict[kmer])
        else:
            score = 0.0
            kmer2 = revcompl(kmer)
            score = float(kmerdict[kmer2])
        kmer_scores_list.append(score)
    return max(kmer_scores_list)


def max_kmer_score_pos(kmerdict, seq):
    """
    Given a dictionary of k-mer scores and a sequence, compute
    the sum scores around  the maximum score

    :param kmerdict:
    :param seq:
    :return:
    """
    k_mers = find_kmers(seq, 8)
    tot_score = []
    for kmer in k_mers:
        if kmer in kmerdict:
            score = float(kmerdict[kmer])
        else:
            score = 0.0
            #kmer2 = revcompl(kmer)
            #score = float(kmerdict[kmer2])
        tot_score.append(score)
    max_pos = tot_score.index(max(tot_score))

    return sum(tot_score[max_pos-4:max_pos+4])


def find_kmers(string, kmer_size):
    """
    Given a  sequence string, extract all k-mers of length kmer_size

    :param string:
    :param kmer_size:
    :return:
    """
    kmers = []
    for i in range(0, len(string)-kmer_size+1):
        kmers.append(string[i:i+kmer_size])
    return kmers


def getkey(item):
    """
    :param item:
    """

    return item[1]
