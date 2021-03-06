import random
import numpy as N
import math
import os.path
import logging
import ext_clustering as ext
import scipy.io as NIO
import Motif

import warnings
warnings.simplefilter("ignore", N.ComplexWarning)

max_iterations = 200
max_repetitions = 10
eps = 0.01
m = 2.0


def getClusters(clustering):
    """
    Transform from one clustering representation to another
    Input:
    - clustering: list l of n elements where l[i] equals the cluster where i belongs
    Output:
    - list with the clusters
    """
    k = max(clustering) + 1
    clusters = [[] for i in range(0, k)]

    i = 0
    for e in clustering:
        clusters[e].append(i)
        i += 1

    return clusters


NO_KERNELIZE = 0
NO_NEGATIVE_EIGENVALUES = 1
QUADRATIC_EIGENVALUES = 2


def kernelize(kernel, type_kernel=NO_NEGATIVE_EIGENVALUES):
    """
    All negative eigenvalues are considered to be 0.
    Input:
    - kernel
    Output:
    - real kernel: semidefinite positive
    """
    import scipy.linalg.decomp as decomp
    import scipy as SCI

    if type_kernel == NO_KERNELIZE:
        return kernel

    w, weigen = decomp.eig(kernel)

    if type_kernel == NO_NEGATIVE_EIGENVALUES:
        w_no_neg = []
        for i in w:
            if i.real > 0:
                w_no_neg.append(i)
            else:
                w_no_neg.append(0)

        k_ = N.dot(N.dot(weigen, SCI.diag(w_no_neg)), N.transpose(weigen))
        print k_
        return k_.astype('f')

    if type_kernel == QUADRATIC_EIGENVALUES:
        return N.dot(kernel, kernel)

    return None


def _readNames(finput, discardFirstLine=True):
    """
    It reads names from first column in kernel matrix.
    """
    if type(finput) == file:
        hinput = finput
    else:
        hinput = open(finput, "r")

    if discardFirstLine:
        hinput.readline()  # Reading first line, dimensions

    names = []

    for line in hinput:
        name = line.split()[0]
        # if name in names: return None	# Duplicated name
        names.append(name)

    return names


def _clusterIndexesToNames(names, clusters, discardFirstLine=True):
    """
    Translates index in kernel matrix to object names.
    """
    list1 = []

    for c in clusters:
        list1.append([])

        for o in c:
            list1[-1].append(names[o])

    return list1


def _writeClusters(clusters, output):
    """
    Writes a cluster set to output
    """
    if type(output) == file:
        ho = output
    else:
        ho = open(output, "w")

    for c in clusters:
        for o in c[0:-1]:
            ho.write(o + "\t")

        if len(c) > 1:
            ho.write(c[-1] + "\n")


def normalize01(matrix):
    max_ = max(matrix.flat)
    min_ = min(matrix.flat)
    if min_ < 0: print "Normalization: Negative entry in matrix"

    matrix = N.subtract(matrix, min_)
    matrix = N.divide(matrix, (max_ - min_))

    return matrix


# SCRIPTS AND MAIN

def usage(argv0):
    """
    Main program usage
    """
    print argv0, "kcmeans sim_matrix:file k:int |"


#	print argv0, "fcmeans sim_matrix:file k:int m:float|"
#	print argv0, "kfcmeans sim_matrix:file k:int m:float|"
#	print argv0, "defuzzify sim_matrix:file output:file"



def script_kcmeans(fileIn, k, fileOut="kcmeans", distance_to_similarity=False, apply_kernel=True, normalize=True):
    """
    Main for kcmeans
    """

    logging.info("script_kcmeans %s %s" % (fileIn, k))
    # loadtxt(fname, dtype=<type 'float'>, comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
    # Load data from a text file.
    # similarities = NIO.read_array(fileIn, separator = '\t', columns = (1,-1), lines = (1,-1) )
    names = _readNames(fileIn)
    rows = len(names) + 1
    cols = ()
    for i in range(1, rows):
        cols = cols + (i,)

    similarities = N.loadtxt(fileIn, delimiter='\t', usecols=cols, skiprows=1)
    names = _readNames(fileIn)
    # print names
    # print similarities
    if normalize: similarities = normalize01(similarities)

    if distance_to_similarity:
        sims = 1 - similarities
    else:
        sims = similarities

    if apply_kernel: sims = kernelize(sims)

    clustering = ext.cmeans_from_similarities(sims, k, max_iterations, max_repetitions)
    clusters = getClusters(clustering)
    namedClusters = _clusterIndexesToNames(names, clusters)
    # print "namedClusters",namedClusters
    _writeClusters(namedClusters, fileOut)

    return 0


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        usage(sys.argv[0])
        sys.exit(-1)

    exec ("out = script_" + sys.argv[1] + "(sys.argv[2:])")

    if out == -1:
        usage(sys.argv[0])
        sys.exit(-1)
