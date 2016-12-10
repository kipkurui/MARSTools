"""
Assess_motif is a python module for assessing the motif Performance using
ChIP-seq.
Requires:
    -> A motif file in meme format
    -> ChIP-seq file (for now) in tab delimited format Chr score sequence
    -> A motif scoring framework to use:
        -gomeroccupancyscore
        -sumoccupancyscore
        -maxoccuopancyscore
        -energyscore
        -sumlogoddsscore
        -maxlogoddsscore
        -amaoccupancyscore
        - A few k-mer test functions
    -> Output file
 Usage:
     python Assess_by_score.py <tf> <scoring_function> <user_motif> <chip_seq_list> <results_folder_path>
"""
import os
import sys
from math import exp

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
from sklearn import metrics

from MARS_Suite.utils import rotate_image
#from back_ups.libkmersvm import *

#################################################################################
# TODO: Get all the functions working well and in generic forms and then
# separate out ChIP-seq and PBM main programs
##################################################################################

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get_motif_from_meme(meme, motif="MOTIF"):
    """
    Extract a motif from meme file given a unique motif
    name and create dictionary for sequence scoring

    Default motif name is keyword MOTIF for single motif files.
    """
    name = ""
    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    flag = 0
    check = 0
    with open(meme, "r") as f1:
        for line in f1:
            if line.startswith('MOTIF'):
                if line.split(" ")[1] == motif:
                    # if str(motif) in line:
                    name = line.split(" ")[1]
                    flag += 1
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                # print line
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    check += 1
        return areapwm, name


def get_motif(meme, motif="MOTIF"):
    """
    Extract a motif from meme file given a unique motif
    name and create dictionary for sequence scoring

    Default motif name is keyword MOTIF for single motif files. 
    """

    pwm_dictionary = {}
    pwm_dictionary["A"] = []
    pwm_dictionary["C"] = []
    pwm_dictionary["G"] = []
    pwm_dictionary["T"] = []
    flag = 0
    check = 0
    with open(meme, "r") as f1:
        for line in f1:
            if str(motif) in line:
                flag += 1
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    pwm_dictionary["A"].append(float(words[0]))
                    pwm_dictionary["C"].append(float(words[1]))
                    pwm_dictionary["G"].append(float(words[2]))
                    pwm_dictionary["T"].append(float(words[3]))
                    check += 1
        return pwm_dictionary, int(w)


def rc_pwm(area_pwm, pwm_len):
    """
    Takes as input the forward pwm and returns a reverse
    complement of the motif
    """

    rcareapwm = {}
    rcareapwm["A"] = []
    rcareapwm["C"] = []
    rcareapwm["G"] = []
    rcareapwm["T"] = []
    for i in range(pwm_len):
        rcareapwm["A"].append(area_pwm["T"][pwm_len - i - 1])
        rcareapwm["C"].append(area_pwm["G"][pwm_len - i - 1])
        rcareapwm["G"].append(area_pwm["C"][pwm_len - i - 1])
        rcareapwm["T"].append(area_pwm["A"][pwm_len - i - 1])
    return rcareapwm


###############################################################################
# Scoring functions
##############################################################################


def gomeroccupancyscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the gomer score
    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        gomer_occupancy = 1
        area_pwm_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(pwm_length - 1, 1, -1):
            prod_gomer = 1
            prod_gomer_rc = 1
            for j in range(pwm_length):
                if j <= i:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                elif (j + i) > len(seq) - 1:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                else:
                    # print "got to else"
                    s = seq[j + i]
                    prod_gomer *= pwm_dictionary[s][j]
                    prod_gomer_rc *= area_pwm_rc[s][j]
            gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
        for i in range(len(seq) - 1):
            prod_gomer = 1
            prod_gomer_rc = 1
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq) - 1:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                else:
                    prod_gomer *= pwm_dictionary[seq[j + i]][j]
                    prod_gomer_rc *= area_pwm_rc[seq[j + i]][j]
            gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
        gomer_occupancy = 1 - gomer_occupancy

        return gomer_occupancy


def sumoccupancyscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum occupancy score.
    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        sum_occupancy = 0
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            occupancy = 1
            occupancy_rc = 1
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    occupancy *= 0.25
                    occupancy_rc *= 0.25
                elif seq[j + i] not in ["A", "C", "G", "T"]:
                    occupancy *= 0.25
                    occupancy_rc *= 0.25
                else:
                    occupancy *= pwm_dictionary[seq[j + i]][j]
                    occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
            sum_occupancy += occupancy + occupancy_rc
        return sum_occupancy / 2


def sumlogoddsscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum of the log odds scores.
    
    This is the scoring approach that is used by MEME Suite
    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        sum_log_odds = 1
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            log_odds = 0
            log_odds_rc = 0
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    log_odds += 0.0
                    log_odds_rc += 0.0
                elif seq[j + i] not in ["A", "C", "G", "T"]:
                    log_odds += 0.0
                    log_odds_rc += 0.0
                else:
                    q = pwm_dictionary[seq[j + i]][j]
                    q_rc = pwm_dictionary_rc[seq[j + i]][j]
                    if q == 0 or (q_rc == 0):
                        q = 0.000000000000000000000000000001  
                        # make this as close to zero as possible
                        q_rc = 0.000000000000000000000000000001
                    else:
                        q = pwm_dictionary[seq[j + i]][j]
                        q_rc = pwm_dictionary_rc[seq[j + i]][j]
                    log_odds += (np.log(q / 0.25) / np.log(2)) * 100
                    log_odds_rc += (np.log(q_rc / 0.25) / np.log(2)) * 100
            sum_log_odds += log_odds + log_odds_rc
        return sum_log_odds / 2


def maxlogoddsscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the max of the log odds scores.

    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        log_odds_list = []
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            log_odds_score = 0
            log_odds_score_rc = 0
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    log_odds_score += 0.0
                    log_odds_score_rc += 0.0
                elif seq[j + i] not in ["A", "C", "G", "T"]:
                    log_odds_score += 0.0
                    log_odds_score_rc += 0.0
                else:
                    q = pwm_dictionary[seq[j + i]][j]
                    q_rc = pwm_dictionary_rc[seq[j + i]][j]
                    if q == 0 or q_rc == 0:
                        q = 0.000000000000000000000000000001  
                        # make this as close to zero as possible
                        q_rc = 0.000000000000000000000000000001
                    else:
                        q = pwm_dictionary[seq[j + i]][j]
                        q_rc = pwm_dictionary_rc[seq[j + i]][j]
                log_odds_score += (np.log(q / 0.25) / np.log(2)) * 100
                log_odds_score_rc += (np.log(q_rc / 0.25) / np.log(2)) * 100
            log_odds_list.append(log_odds_score)
            # FIXME: There was an error here in which we did not include 
            #the reverse complement in the computation
            log_odds_list.append(log_odds_score_rc)
        max_log_odds = max(log_odds_list)
        return max_log_odds


def maxoccupancyscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the sum occupancy score.
    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        occupancy_list = []
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            occupancy = 1
            occupancy_rc = 1
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    occupancy *= 0.25
                    occupancy_rc *= 0.25
                elif seq[j + i] not in ["A", "C", "G", "T"]:
                    occupancy *= 0.25
                    occupancy_rc *= 0.25
                else:
                    occupancy *= pwm_dictionary[seq[j + i]][j]
                    occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
            occupancy_list.append(occupancy)
            occupancy_list.append(occupancy_rc)
        max_occupancy = max(occupancy_list)
        return max_occupancy


def amaoccupancyscore(pwm_dictionary, seq):
    """
    Score sequences using AMA scoring uses average occupancy scores

    """
    if "N" in seq:
        return 0
    else:
        # pwm_length = len(pwm_dictionary)
        pwm_length = len(pwm_dictionary["A"])
        occupancy_list = []
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            occupancy = 1
            occupancy_rc = 1
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    occupancy *= 0.25
                    occupancy_rc *= 0.25
                else:
                    occupancy *= pwm_dictionary[seq[j + i]][j]
                    occupancy_rc *= pwm_dictionary_rc[seq[j + i]][j]
            occupancy_list.append(occupancy + occupancy_rc)
        ama_occupancy = sum(occupancy_list) / len(occupancy_list)
        return ama_occupancy


###############################################################################
# TODO convert the scoring function into class object to reduce
##############################################################################


def energyscore(pwm_dictionary, seq):
    """
    Score sequences using the beeml energy scoring approach.

    Borrowed greatly from the work of Zhao and Stormo

    P(Si)=1/(1+e^Ei-u)

    Ei=sumsum(Si(b,k)e(b,k))

    Previous approaches seem to be using the the minimum sum of the
    energy contribution of each of the bases of a specific region.

    This is currently showing some promise but further testing is
    needed to ensure that I have a robust algorithm.
    """
    if "N" in seq:
        return 0
    else:
        pwm_length = len(pwm_dictionary["A"])
        energy_list = []
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            energy = 0
            energy_rc = 0
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    energy += 0.25
                    energy_rc += 0.25
                else:
                    energy += pwm_dictionary[seq[j + i]][j]
                    energy_rc += pwm_dictionary_rc[seq[j + i]][j]

                energy_list.append(1 / (1 + (exp(energy))))
                energy_list.append(1 / (1 + (exp(energy_rc))))
        energy_score = min(energy_list)
        return energy_score


######################################################################
#      Assessment metrics
######################################################################
# Convert all these assessment metrics into a class object and do the same
# for the rest of the  functions to respective categories of classes.
#####################################################################

# TODO: Replace auc with this when I start using the pandas dataframe


def compute_auc(predicted, cutoff, label):
    """
    Compute Area under ROC curve given a prediction file formatted as floows:
        seq_name                \t Score           \t Classification \n
        chr2:43019807-43019857	\t 4.251985e-06	\t 1 \n
        . \n
        . \n
        . \n
        chr2:43619807-43619857	\t 4.251985e-08	\t 0 \n
    log: Changed to be able to compute AUC for any scores
    """
    y = np.concatenate((np.ones(cutoff), np.zeros(cutoff)), axis=0)

    fpr, tpr, thresholds = metrics.roc_curve(y, predicted[:cutoff*2], pos_label=label)

    auc = metrics.auc(fpr, tpr)

    return auc


def compute_auc_old(predicted, cutoff, label):
    """
    Compute Area under ROC curve given a prediction file formatted as floows:
        seq_name                \t Score           \t Classification \n
        chr2:43019807-43019857	\t 4.251985e-06	\t 1 \n
        . \n
        . \n
        . \n
        chr2:43619807-43619857	\t 4.251985e-08	\t 0 \n
    log: Changed to be able to compute AUC for any scores
    """
    y = np.concatenate((np.ones(cutoff), np.zeros(cutoff)), axis=0)

    pr = np.array(predicted)

    fpr, tpr, thresholds = metrics.roc_curve(y, pr, pos_label=label)
    auc = metrics.auc(fpr, tpr)
    return auc


def compute_mncp(predicted, cutoff, label):
    """
    This is the MNCP computation adopted from Clarke 2003

    MNCP is a rank based metric similar to AUC but 
    its a plot of TP and all positives
    hence considered to be less affected by false positives.

    MNCP is the mean normalized

    Functions adapted from GimmeMotifs tools by Simon Van Heeringen
    """
    from numpy import mean, array, hstack
    if label == 1:
        fg_vals = predicted[:cutoff]
        bg_vals = predicted[cutoff:]
    else:
        fg_vals = predicted[cutoff:]
        bg_vals = predicted[:cutoff]
    fg_len = len(fg_vals)
    total_len = len(fg_vals) + len(bg_vals)

    if type(fg_vals) != type(array([])):
        fg_vals = array(fg_vals)
    if type(bg_vals) != type(array([])):
        bg_vals = array(bg_vals)
    # Rank the data
    fg_rank = stats.rankdata(fg_vals)

    # combine foreground and background data and get the ranks
    total_rank = stats.rankdata(hstack((fg_vals, bg_vals)))
    slopes = []
    for i in range(len(fg_vals)):
        slope = ((fg_len - fg_rank[i] + 1) / fg_len) / ((total_len - total_rank[i] + 1) / total_len)
        slopes.append(slope)
    mncp = mean(slopes)
    return mncp


def compute_spearman(observed, predicted, cutoff, label=0):
    """
    Compute Spearman correlation in the positive set of ChIP-seq data.

    Given ChIp-seq data, we need to compute the correlation of the predicted
    score and the intensity score
    :param label:
    """

    speaman = stats.spearmanr(observed[:cutoff], predicted[:cutoff])[0]

    if label == 0:
        speaman *= -1

    return speaman


def compute_pearson(observed, predicted, cutoff, label=1):
    """
    Compute pearson correlation in the positive set of data.

    Given data, we need to compute the correlation of the predicted
    score and the existing score
    """
    observed = [float(i) for i in observed]
    predicted = [float(i) for i in predicted]
    observed = [float(i) - np.mean(observed) for i in observed]
    predicted = [float(i) - np.mean(predicted) for i in predicted]

    pearson = stats.pearsonr(observed[:cutoff], predicted[:cutoff])[0]
    if label == 0:
        pearson *= -1
    return pearson


############################################################
# Assess motifs using ChIP-seq data
############################################################


def score_chipseq(chip_seq, score_function, user_motif_details):
    """
    Give the ChIP-seq file, this function scores the sequences
    using the specified scoring function.
    """

    # if user_motif_details:
    area_pwm = user_motif_details[0]
    # pwm_length = len(area_pwm['A'])

    test_file = pd.read_table(chip_seq, header=None)
    if len(test_file.columns) == 3:
        seq_col = 2
        chip_score = test_file[1]
    else:
        chip_score = np.zeros(len(test_file))
        seq_col = 1

    seq_score = test_file[seq_col].apply(lambda seq: score_function(area_pwm, seq))

    return chip_score, seq_score


#############################################################################
# Run the complete assess program as a function
#############################################################################


def run_assess(score_function, summary_output, raw_output, user_motif_details, chip_seq_list):
    with open(summary_output, "a") as out:
        auc = []
        spearman = []
        mncp = []
        pearson = []
        score = eval(score_function)
        if score_function == "energyscore":  # Lower values are better, hence the reversing in Energy
            label = 0
        else:
            label = 1
        motif = user_motif_details[1].strip()

        # Extract the cell line and lab name from the ChIP-seq file
        with open(raw_output, 'a') as raw_out:
            for raw_chip_data in chip_seq_list:
                cell_lab = raw_chip_data.split('/')[-1].split('.')[0]

                chip_score = score_chipseq(raw_chip_data, score, user_motif_details)
                cut_off = len(chip_score[1]) / 2  # use a flexible cut-off dictated by the sze of the input file

                au = compute_auc(chip_score[1], cut_off, label)
                auc += [au]
                mn = compute_mncp(chip_score[1], cut_off, label)
                mncp += [mn]
                sp = compute_spearman(chip_score[0], chip_score[1], cut_off, label)
                spearman.append(sp)
                pe = compute_pearson(chip_score[0], chip_score[1], cut_off, label)
                pearson.append(pe)

                write_out_data = "%s\t%s\t%f\t%f\t%f\t%f\n" % (cell_lab, motif, au, mn, pe, sp)
                raw_out.write(write_out_data)

            # Write out the raw data
            write_out_data = "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % ('Average', motif, np.mean(auc), np.mean(mncp),
                                                                   np.mean(pearson), np.mean(spearman))
            raw_out.write(write_out_data)

        # Write out teh summary data
        write_out_data = "%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (motif, np.mean(auc), np.mean(mncp), np.mean(pearson),
                                                           np.mean(spearman))
        out.write(write_out_data)


score_extensions = {"gomeroccupancyscore": 'gomer', "energyscore": 'energy', "amaoccupancyscore": 'ama',
                    "maxoccupancyscore": 'maxoc', "sumlogoddsscore": 'sumlog', "maxlogoddsscore": 'maxlog',
                    "sumoccupancyscore": 'sumoc', "energy_score_kmer": 'energymer', "kmersvm_score_seq": 'kmersvm',
                    "max_score_kmer": 'max_kmer', "max_score_kmer_pos": 'max_kmer_pos'}


def get_tf_names(user_motif):
    tf_names = []
    with open(user_motif) as motif_input:
        for line in motif_input:
            if line.startswith('MOTIF'):
                mot_name = line.split(" ")[1]
                tf_names.append(mot_name)
    return tf_names


def run_all(tf, scoring_function, user_motif, chip_seq_list, results_folder_path):
    import random
    if len(chip_seq_list) > 10:
        random.seed(10)
        chip_seq_list = random.sample(chip_seq_list, 10)

    score_option = scoring_function
    summary_output_file = "%s/%s.%s" % (results_folder_path, tf.lower(), score_extensions[scoring_function])
    raw_output_file = "%s/%s_raw.%s" % (results_folder_path, tf.lower(), score_extensions[scoring_function])
    pr = "%s\t%s\t%s\t%s\t%s\n" % ("Motif", "AUC", "MNCP", "Pearson", "Spearman")

    file_header = "%s\t%s\t%s\t%s\t%s\t%s\n" % ("Cell_lab", "Motif", "AUC", "MNCP", "Pearson", "Spearman")

    with open(summary_output_file, "w") as write_summary:
        with open(raw_output_file, "w") as write_raw:
            write_summary.write(pr)
            write_raw.write(file_header)
    if type(user_motif) is list:

        # A crude way of determining the assessment i'm running
        
        tf_names = user_motif
        for mot_name in tf_names:
            global g_kmer2id
            g_kmer2id = get_g_kmer2id(8)
            kmer_name = mot_name.split("/")[-1]
            if score_option == "kmersvm_score_seq":
                kmer_svmw_dict = get_kmer_dict_rev(mot_name, kmer_name)
                user_motif_details = get_svm_list(kmer_svmw_dict[0], kmer_name)

            else:
                user_motif_details = get_kmer_dict_rev(mot_name, kmer_name)

            run_assess(score_option, summary_output_file, raw_output_file, user_motif_details, chip_seq_list)
        
    else:
        tf_names = get_tf_names(user_motif)
        for mot_name in tf_names:
            user_motif_details = get_motif_from_meme(user_motif, mot_name)
            run_assess(score_option, summary_output_file, raw_output_file, user_motif_details, chip_seq_list)

        # call the plotting function
        plot_info(tf, scoring_function, results_folder_path)
###############################################################################
#
#    kmer scoring stuff
#
#############################################################################


def get_kmer_dict_rev(kmerscore, kmer_name):
    test = pd.read_table(kmerscore, index_col="8-mer", usecols=["8-mer", "E-score"])
    test.fillna(0, inplace=True)
    test2 = pd.read_table(kmerscore, index_col="8-mer.1", usecols=["8-mer.1", "E-score"])
    test2.index.name = "8-mer"
    test2.fillna(0, inplace=True)
    combined = test.append(test2)
    combined_dict = combined.to_dict()["E-score"]

    return combined_dict, kmer_name


# def get_g_kmer2id(k):
#     g_kmer2id = {}
#     kmers = generate_kmers(k)
#     rcmap = generate_rcmap_table(k, kmers)
#     for i in xrange(len(kmers)):
#         g_kmer2id[kmers[i]] = rcmap[i]
#     for i in xrange(len(kmers)):
#         g_kmer2id[kmers[i]] = rcmap[i]
#
#     return g_kmer2id


def get_svm_list(kmer_svmw_dict, name):
    # svmw_list = []
    svmw = [0]*(2**(2*8))

    for kmer in kmer_svmw_dict.keys():
        svmw[g_kmer2id[kmer]] = kmer_svmw_dict[kmer]

    svmw_list = np.array(svmw, np.double)

    return svmw_list, name


def kmersvm_score_seq(kmer_svmw_dict, seq):
    """calculate SVM score of given sequence using single set of svm weights

    Arguments:
    s -- string, DNA sequence
    kmer_svmw_dict -- A dictionary of kmers
    kmerlen -- integer, length of k-mer of SVM weight

    Return:
    SVM score
    """
    kmerlen = 8
    kmer2id = g_kmer2id
    x = [0]*(2**(2*kmerlen))
    for j in xrange(len(seq)-kmerlen+1):
        x[kmer2id[seq[j:j+kmerlen]]] += 1

    x = np.array(x, np.double)
    score_norm = np.dot(kmer_svmw_dict, x)/np.sqrt(np.sum(x**2))

    return score_norm

revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in x][::-1])


def score_kmer(kmer, kmerdict, revcompl):
    score = 0
    if kmer in kmerdict:
        score = float(kmerdict[kmer])
    else:
        kmer2 = revcompl(kmer)
        score = float(kmerdict[kmer2])
    return score


def energy_score_kmer(kmerdict, seq):
    k_mers = find_kmers(seq, 8)
    tot_score = 0
    for kmer in k_mers:
        if kmer in kmerdict:
            score = float(kmerdict[kmer])
        else:
            kmer2 = revcompl(kmer)
            score = float(kmerdict[kmer2])
        tot_score += score
    return tot_score


def max_score_kmer(kmerdict, seq):
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
    return max(tot_score)


def max_score_kmer_pos(kmerdict, seq):
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
    kmers = []
    for i in range(0, len(string)-kmer_size+1):
        kmers.append(string[i:i+kmer_size])
    return kmers


def getkey(item):
    return item[1]


# def get_kmer_dict_rev(kmerscore, kmer_name):
#
#     test = pd.read_table(kmerscore, index_col="8-mer", usecols=["8-mer", "E-score"])
#     test = test[test["E-score"] >0.35]
#     scoredict = test.to_dict()["E-score"]
#     # with open(kmerscore) as kmers:
#     #     print kmerscore
#     #     for line in kmers:
#     #         ke, rems, val = line.split("\t")
#     #
#     #         scoredict[ke] = val
#     return scoredict, kmer_name


def get_kmer_dict(kmerscore, kmer_name):
    """
    Convert the forward and reverse kmers into a dictionary
    for a quick look-up and scoring of sequences
    """
    test = pd.read_table(kmerscore, index_col="8-mer", usecols=["8-mer", "E-score"])
    test.fillna(0, inplace=True)
    test2 = pd.read_table(kmerscore, index_col="8-mer.1", usecols=["8-mer.1", "E-score"])
    test2.index.name = "8-mer"
    test2.fillna(0, inplace=True)
    combined = test.append(test2)
    combined_dict = combined.to_dict()["E-score"]

    return combined_dict, kmer_name
    
##############################################################################
# Plotting
###############################################################################


def plot_info(tf, score_method, results_folder):
    files_path = '%s/%s' % (results_folder, tf)
    score_ext = score_extensions[score_method]

    plot_histogram_assess(files_path + "." + score_ext,
                          files_path + '_assess')
    plot_histogram_assess(files_path + "." + score_ext,
                          files_path + '_assess.eps')
    plot_raw_assess(files_path + "_raw." + score_ext, files_path + '_assess_raw.png', 'AUC')
    plot_raw_assess(files_path + "_raw." + score_ext, files_path + '_assess_raw.eps', 'AUC')
    rotate_image(files_path + '_assess.png', files_path + '_assess_rot.png')


def plot_raw_assess(raw_data, figure_output, stat):
    """
    This function allows the raw data to be plotted  in the form of a clustermap

    This way, information about how each motif scores in different cell lines is
    obtained
    """
    sns.set_style("white")
    raw_max = pd.read_table(raw_data)
    raw_max = raw_max.drop_duplicates()

    raw_edit = raw_max.pivot('Motif', 'Cell_lab', stat)
    raw_edit.sort(columns="Average", axis=0, ascending=False, inplace=True)
    cg = sns.clustermap(raw_edit, method='single', metric="euclidean", z_score=None,
                        annot=True, row_cluster=False, col_cluster=True, linewidths=.15)
    # to rotate the y-axis labels correctly
    test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    f = plt.gcf()
    f.savefig(figure_output, bbox_inches='tight')


def plot_histogram_assess(assess_input, figure_output):
    """
    This function plots all the assess motifs output on a single plot and
    sorts the data based on a statistical function used.

    The function could be modified for flexibility in that it can be used
    to generate both plotly and seaborn figures. Also, a choice of the
    statistical tool used for sorting the data can be added through the
    parameters.

    Borrow something from this later on:
    http://geophysik.uni-muenchen.de/~krischer/___bla___/S
    eaborn.docset/Contents/Resources/Documents/stanford.edu/
    _mwaskom/software/seaborn/tutorial/axis_grids.html
    """
    # print "We got to plot histogram"
    sns.set_style("white")
    raw_auc = pd.read_table(assess_input, index_col="Motif")
    raw_auc = raw_auc.drop_duplicates()
    # df = df.T.drop_duplicates().T
    raw_auc = raw_auc.sort(columns="MNCP", axis=0, ascending=False)
    labels = raw_auc.index
    x = 10
    if len(labels) > 50:
        x = 15
    elif len(labels) < 10:
        x = 5
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(x, 10), sharex=True)
    a = sns.barplot(x=labels, y=raw_auc["AUC"],
                    palette='colorblind', x_order=labels, ax=ax1)
    b = sns.barplot(x=labels, y=raw_auc["MNCP"],
                    palette="colorblind", x_order=labels, ax=ax2)
    c = sns.barplot(x=labels, y=raw_auc["Spearman"],
                    palette="colorblind", x_order=labels, ax=ax3)
    d = sns.barplot(x=labels, y=raw_auc["Pearson"],
                    palette="colorblind", x_order=labels, ax=ax4)
    d.set_xticklabels(labels, rotation=90)

    sns.despine()
    f.savefig(figure_output + ".eps", bbox_inches='tight')
    f.savefig(figure_output + ".png", bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print __doc__
        sys.exit(1)
    tf_par = sys.argv[1]
    scoring_fun = sys.argv[2]
    user_mot = sys.argv[3]
    chip_list = sys.argv[4]
    results_path = sys.argv[5]

    run_all(tf_par, scoring_fun, user_mot, chip_list, results_path)



