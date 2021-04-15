from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from .Assess_by_score import np, sys, SCORE_EXTENSION, compute_mncp, \
    compute_pearson, compute_spearman,\
    get_motif_from_meme, pd, compute_auc


###########################################################################
# Assess motifs using PBM data
############################################################################


def get_pbm_possitives(pbm):
    """
    Normalize Pbm data to extract positive and negative probes
    for AUC computation.

    Borrowed greatly from Cheng 2007

    The PBM intensity score are transformed to positive sequences by adding
    the minimum score+1 to each score then log transformed.

    The positive probes are generally, 4*MAD(mean absolute deviation) above
    the median of the normalized scores.

    The rest are considered as the negative probes but may revise this.
    """

    ##TODO: Include an option for optional printout of the output
    intensity = []
    positive = []
    with open(pbm) as probes:
        for line in probes:
            intensity.append(float(line.split()[0]))
            #sequences.append(line.split()[1])
        intensity = np.array(intensity)
        minimum = np.min(intensity) - 1
        intensity -= (minimum)
        transformed = np.log(intensity)
    median_av = np.median(transformed) + (0.6745 * 4)

    return np.sum(transformed > median_av)


# def score_pbm(pbm, score_function, user_motif_details):
#     """
#     Given the PBM probe file, this function scores the sequences
#     using the specified scoring function
#
#     :param pbm:
#     :param score_function:
#     :param user_motif_details:
#     :return:
#     """
#     pbm_score = []
#     intensity = []
#     pos = get_pbm_possitives(pbm)
#     if pos < 50:
#         pos = 50
#     os.system("head -%i %s >/tmp/pos.txt" % (pos, pbm))
#     os.system("tail -%i %s >>/tmp/pos.txt" % (pos, pbm))
#     with open("/tmp/pos.txt") as debru:
#         for line in debru:
#             details = line.split()
#             seq = details[1]
#             if user_motif_details:
#                 area_pwm = user_motif_details[0]
#             pwm_length = len(area_pwm['A'])
#             seq_score = score_function(area_pwm, pwm_length, seq[:36])
#             pbm_score.append(seq_score)
#             intensity.append(details[0])
#     os.system("rm /tmp/pos.txt")
#     return intensity, pbm_score

def score_pbm(pbm, score_function, user_motif_details):
    """
    Given the PBM probe file, this function scores the sequences
    using the specified scoring function

    :param pbm:
    :param score_function:
    :param user_motif_details:
    :return:

    This function now does what was being done by two
    """
    pbm_score = []
    intensity = []

    pbms = pd.read_table(pbm, header=None)

    transformed = np.log(pbms[0] - (np.min(pbms[0])-1))

    pos = np.sum(transformed > (np.median(transformed) + (0.6745 * 4)))

    if pos < 50:
        pos = 50

    pbm_use = pbms.head(pos).append(pbms.head(pos))

    seq_col = 1
    area_pwm = user_motif_details
    pbm_score = pbm_use[seq_col].apply(lambda seq: score_function(area_pwm, seq[:36]))

    intensity = pbm_use[0]

    return intensity, pbm_score
##############################################################################
# Run the complete assess program as a function
#################################################################################


def run_assess_pbm(score_function, summary_output, raw_output, user_motif_details, pbm_list):
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

        with open(raw_output, 'a') as raw_out:
            for raw_pbm_data in pbm_list:
                if get_pbm_possitives(raw_pbm_data) == 0:
                    continue
                else:

                    cell_lab = raw_pbm_data.split('/')[-1].split('_deBruijn')[0]  # Specific to PBM
                    pbm_score = score_pbm(raw_pbm_data, score, user_motif_details)  # Specific to PBM
                    cut_off = len(pbm_score[1]) // 2  # use a flexible cut-off dictated by the sze of teh input file

                    au = compute_auc(pbm_score[1], cut_off, label)
                    auc += [au]
                    mn = compute_mncp(pbm_score[1], cut_off, label)
                    mncp += [mn]
                    sp = compute_spearman(pbm_score[0], pbm_score[1], cut_off, label)
                    spearman.append(sp)
                    pe = compute_pearson(pbm_score[0], pbm_score[1], cut_off, label)
                    pearson.append(pe)

                    write_out_data = "%s\t%s\t%f\t%f\t%f\t%f\n" % (cell_lab, motif, au, mn, pe, sp)
                    raw_out.write(write_out_data)

                # Write out the raw data
            write_out_data = "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % ('Average', motif, np.mean(auc), np.mean(mncp),
                                                                       np.mean(pearson), np.mean(spearman))
            raw_out.write(write_out_data)

        # Write out the summary data
        write_out_data = "%s\t%.3f\t%.3f\t%.3f\t%.3f\n" % (motif, np.mean(auc), np.mean(mncp), np.mean(pearson),
                                                           np.mean(spearman))
        out.write(write_out_data)


def get_tf_names(user_motif):
    tf_names = []
    with open(user_motif) as motif_input:
        for line in motif_input:
            if line.startswith('MOTIF'):
                mot_name = line.split(" ")[1]
                tf_names.append(mot_name)
    return tf_names


def run_all_pbm(tf, scoring_function, user_motif, pbm_list, results_folder_path):
    import random
    if len(pbm_list) > 10:
        random.seed(10)
        pbm_list = random.sample(pbm_list, 10)

    score_option = scoring_function
    summary_output_file = "%s/%s.%s" % (results_folder_path, tf.lower(), score_extensions[scoring_function])
    raw_output_file = "%s/%s_raw.%s" % (results_folder_path, tf.lower(), score_extensions[scoring_function])
    pr = "%s\t%s\t%s\t%s\t%s\n" % ("Motif", "AUC", "MNCP", "Pearson", "Spearman")

    file_header = "%s\t%s\t%s\t%s\t%s\t%s\n" % ("Cell_lab", "Motif", "AUC", "MNCP", "Pearson", "Spearman")

    with open(summary_output_file, "w") as write_summary:
        with open(raw_output_file, "w") as write_raw:
            write_summary.write(pr)
            write_raw.write(file_header)

    tf_names = get_tf_names(user_motif)

    for mot_name in tf_names:
        user_motif_details = get_motif_from_meme(user_motif, mot_name)
        run_assess_pbm(score_option, summary_output_file, raw_output_file,\
             user_motif_details, pbm_list)


# TODO: Add a plotting function independent in each module
if __name__ == '__main__':
    if len(sys.argv) < 6:
        print(__doc__)
        sys.exit(1)
    tf = sys.argv[1]
    scoring_function = sys.argv[2]
    user_motif = sys.argv[3]
    pbm_list = sys.argv[4]
    results_folder_path = sys.argv[5]

    run_all_pbm(tf, scoring_function, user_motif, pbm_list, results_folder_path)
