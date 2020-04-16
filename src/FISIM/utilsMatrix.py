from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div
from .Motif import *
import numpy as np


def computeSimMatrix(motifs, cores, distance=Motif.fisim, verbose=False):
    mat = np.zeros((len(motifs)+1, len(motifs)+1), dtype=float)
    if verbose:
        print("Computando la similitud entre " + str(len(motifs)) + " motivos")
    for i in range(len(motifs)):
        if verbose:
            print("Procesando motivo " + str(i+1))
        for j in range(len(motifs)):
            if len(cores) == 0:
                mat[i, j] = distance(motifs[i], motifs[j])[0]
            else:
                mat[i, j] = old_div((distance(cores[i], motifs[j])[0] + distance(motifs[i], cores[j])[0]), 2)
    return mat


def matrix2File(matrix, motifs, fileOut, ID=False,):
    """
    :param matrix:
    :param motifs:
    :param fileOut:
    :param figout:
    :param ID:
    :return:
    """
    #z= N.zerosmatrix
    f = open(fileOut, 'w')
    rows = []
    old_matrix = matrix
    for i in range(len(motifs)):
        rows.append(motifs[i].ID)
        matrix[i][-1] = old_div(sum(matrix[i]),len(motifs))
    matrix = matrix.T
    for i in range(len(motifs)+1):
        #rows.append(motifs[i].ID)
        matrix[i][-1] = old_div(sum(matrix[i]),len(motifs))
    for i in range(len(motifs)):
        rows.append(motifs[i].ID)
    matrix = matrix.T
    rows.append("Average")
    matrix = matrix.tolist()
    for i in range(len(motifs)+1):
        matrix[i] = [rows[i]] + matrix[i]
    matrix = np.array(matrix)
    matrix = matrix[np.argsort(matrix[:, -1])]
    rows = []
    for i in range(len(matrix)):
        rows.append(matrix[i][0])
    f.write("Motifs")
    for i in range(len(motifs)):
        f.write('\t'+motifs[i].name)
    f.write('\t'+"Average")
    for i in range(len(motifs)):
        #f.write('\t'+motifs[i].ID)
        if ID:
            f.write('\n'+motifs[i].ID)  # primera columna con los IDs
        else:
            f.write('\n'+motifs[i].name)  # primera columna con los nombres
        for j in range(len(motifs)+1):
            f.write('\t' + str(old_matrix[i, j]))
    f.write('\n')  # El \n final


def createSimMatrix(motifs, cores, fileOut, distance=Motif.fisim,ID=False,verbose=False):
    mat = computeSimMatrix(motifs, cores, distance=distance, verbose=verbose)
    matrix2File(mat, motifs, fileOut, ID=ID)