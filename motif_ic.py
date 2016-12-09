#!/usr/bin/python


from __future__ import print_function
import sys
from math import log
import os


# TODO: Give this function a sensible name
def get_motif_summary_tfid(motif_file, out_file, meme_path, reults_folder):

    test = motif_file
    motif_file2 = test
    found = 0
    row = 0
    n_rows = 0
    tn_rows = 0
    entropy = 0
    total_entropy = 0
    motifs = 0
    name = ""
    n = 0
    raw_dict = {}
    #job_no = reults_folder.split("/")[-1]
    with open(out_file, "w") as write_out:

        # with open(raw_file) as raw_in:
        #     for line in raw_in:
        #         raw_dict[line.split()[0]] = line.split()[-1]
        out = "Motif_name\tMotif_IC\tAverage_IC\tMotif_length\tMotif_logo\n"
        write_out.write(out)
        with open(motif_file, "r") as motif_file:
            for line in motif_file:
                words = line.split()
                if found == 0:
                    if line.startswith("MOTIF"):
                        # allow for motifs without an alternative name
                        if len(words) < 3:
                            words.append("")
                        name = (words[1])
                        found = 1
                        motifs += motifs
                        entropy = 0
                        continue
                if found == 1:
                    if line.startswith("letter-probability"):
                        n_rows = int((line.split("w="))[1].split()[0])
                        found = 2
                    continue
                if found == 2:
                    if line == "\n":
                        continue
                    else:
                        check = 0
                    for val in words:
                        if float(val) > 0:
                            check += float(val) * log(float(val))/log(2.0)
                            entropy += float(val) * log(float(val))/log(2.0)
                    row += 1
                    if row >= n_rows:
                        v = 2*n_rows+entropy
                        out = '%s\t%f\t%f\t%i\t<img src="/static/files/temp/%s.png" alt="My image" class="img-responsive"/>\n'\
                              % (name, v, (v/n_rows), n_rows, name)
                        write_out.write(out)
                        found = 0
                        row = 0
                        total_entropy += (v/n_rows)

        mot_list = []
        #os.system("mkdir -p %s/temp" % (reults_folder, tf_name))
        with open(test) as meme_in:
            for line in meme_in:
                if line.startswith("MOTIF"):
                    mot_list.append(line.split()[1])
        fold = "%s" % reults_folder
        for i in mot_list:
            # TODO: Check at this point if data is available in the database, to ensure no repeat discovery
            os.system("%s/ceqlogo -i%s %s -Y -f PNG -h3 -w4 -o %s/%s.png" % (meme_path, i, motif_file2, fold, i))


def get_motif_summary(motif_file, raw_file, tf_name, out_file, meme_path, reults_folder):

    test = motif_file
    motif_file2 = test
    found = 0
    row = 0
    n_rows = 0
    tn_rows = 0
    entropy = 0
    total_entropy = 0
    motifs = 0
    name = ""
    n = 0
    raw_dict = {}
    job_no = reults_folder.split("/")[-1]
    with open(out_file, "w") as write_out:

        with open(raw_file) as raw_in:
            for line in raw_in:
                raw_dict[line.split()[0]] = line.split()[-1]
        out = "Motif_name\tMotif_IC\tAverage_IC\tMotif_length\tMotif_score\tMotif_logo\n"
        write_out.write(out)
        with open(motif_file, "r") as motif_file:
            for line in motif_file:
                words = line.split()
                if found == 0:
                    if line.startswith("MOTIF"):
                        # allow for motifs without an alternative name
                        if len(words) < 3:
                            words.append("")
                        name = (words[1])
                        found = 1
                        motifs += motifs
                        entropy = 0
                        continue
                if found == 1:
                    if line.startswith("letter-probability"):
                        n_rows = int((line.split("w="))[1].split()[0])
                        found = 2
                    continue
                if found == 2:
                    if line == "\n":
                        continue
                    else:
                        check = 0
                    for val in words:
                        if float(val) > 0:
                            check += float(val) * log(float(val))/log(2.0)
                            entropy += float(val) * log(float(val))/log(2.0)
                    row += 1
                    if row >= n_rows:
                        v = 2*n_rows+entropy
                        out = '%s\t%f\t%f\t%i\t%f\t<img src="/static/files/compare/%s/%s/motifs/%s.png" alt="My image" class="img-responsive"/>\n'\
                              % (name, v, (v/n_rows), n_rows, float(raw_dict[name]), job_no, tf_name, name)
                        write_out.write(out)
                        found = 0
                        row = 0
                        total_entropy += (v/n_rows)

        mot_list = []
        os.system("mkdir -p %s/%s/motifs" % (reults_folder, tf_name))
        with open(test) as meme_in:
            for line in meme_in:
                if line.startswith("MOTIF"):
                    mot_list.append(line.split()[1])
        fold = "%s/%s/motifs" % (reults_folder, tf_name)
        for i in mot_list:
            os.system("%s/ceqlogo -i%s %s -Y -f PNG -h3 -w4 -o %s/%s.png" % (meme_path, i, motif_file2, fold, i))


def motif_summary(motif_file, raw_file, out_file):
    found = 0
    row = 0
    n_rows = 0
    entropy = 0
    total_entropy = 0
    motifs = 0
    name = ""
    raw_dict = {}
    with open(out_file, "w") as write_out:

        with open(raw_file) as raw_in:
            for line in raw_in:
                raw_dict[line.split()[0]] = line.split()[1:5]
        out = "Motif_name\tMotif_IC\tAverage_IC\tMotif_length\t%s\t%s\t%s\t%s\n" % \
              (raw_dict["Motif"][0], raw_dict["Motif"][1], raw_dict["Motif"][2], raw_dict["Motif"][3])

        write_out.write(out)
        with open(motif_file, "r") as motif_file:
            for line in motif_file:
                words = line.split()
                if found == 0:
                    if line.startswith("MOTIF"):
                        # allow for motifs without an alternative name
                        if len(words) < 3:
                            words.append("")
                        name = (words[1])
                        found = 1
                        motifs += motifs
                        entropy = 0
                        continue
                if found == 1:
                    if line.startswith("letter-probability"):
                        n_rows = int((line.split("w="))[1].split()[0])
                        found = 2
                    continue
                if found == 2:
                    if line == "\n":
                        continue
                    else:
                        check = 0
                    for val in words:
                        if float(val) > 0:
                            check += float(val) * log(float(val))/log(2.0)
                            entropy += float(val) * log(float(val))/log(2.0)
                    row += 1
                    if row >= n_rows:
                        v = 2*n_rows+entropy
                        out = '%s\t%f\t%f\t%i\t%f\t%f\t%f\t%f\n'\
                              % (name, v, (v/n_rows), n_rows, float(raw_dict[name][0]), float(raw_dict[name][1]),
                                 float(raw_dict[name][2]), float(raw_dict[name][3]))
                        write_out.write(out)
                        #n+= 1
                        #print(n)
                        found = 0
                        row = 0
                        total_entropy += (v/n_rows)