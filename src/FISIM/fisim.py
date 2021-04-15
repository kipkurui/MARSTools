from __future__ import print_function
from __future__ import absolute_import
import utilsMatrix
from Motif import *
import sys
import os

import utils


def helpProgram():
    print("\nCorrect usage:\n\n\tpython fisim.py -f1 <FileIn1> -f2 <FileIn2> [-o <fileOut>] [-v]")
    print("or\n\tpython fisim.py -fileList <FileIn> -o <FileOut> [-core] [-ID]")
    print("\nFirst option: fisim will be computed for the two motifs found in f1 and f2.")
    print("\nSecond option: fisim will be computed pairwise for all the motifs found in fileList. The output will be a similarity matrix.")
    print("\t-core is used when one wants to compute similarities between the cores.")
    print("\t-ID is used when ones want to use the ID instead of the name for motifs identification.")
    print("\npython fisim -h prints this help\n")

f1 = None
f2 = None
fileList = None
out = None
verbose = False
core = False
ID = False
fig_out = None

i = 1
while i < len(sys.argv):
    if sys.argv[i] == "-f1":
        i += 1
        f1 = sys.argv[i]
    elif sys.argv[i] == "-f2":
        i += 1
        f2 = sys.argv[i]
    elif sys.argv[i] == "-fileList":
        i += 1
        fileList = sys.argv[i]
    elif sys.argv[i] == "-o":
        i += 1
        out = sys.argv[i]
    elif sys.argv[i] == "-core":
        core = True
    elif sys.argv[i] == "-ID":
        ID = True
    elif sys.argv[i] == "-v":
        verbose = True
    elif sys.argv[i] == "-h":
        helpProgram()
        sys.exit(-1)
    else:
        print("Unknown flag " + sys.argv[i])
        helpProgram()
        sys.exit(-1)
    i += 1

if f1 is not None and f2 is None:
    print("-f2 not found")
    helpProgram()
    sys.exit(-1)

if f1 is None and f2 is not None:
    print("-f1 not found")
    helpProgram()
    sys.exit(-1)

if fileList is not None and out is None:
    print("Output file not found")
    helpProgram()
    sys.exit(-1)


if f1 is None and f2 is None and fileList is None:
    print("Input files not found")
    helpProgram()
    sys.exit(-1)

if f1 is not None and not os.path.isfile(f1):
    print("Input file f1 " + f1 + " not found")
    sys.exit(-1)

if f2 is not None and not os.path.isfile(f2):
    print("Input file f2 " + f2 + " not found")
    sys.exit(-1)

if fileList is not None and not os.path.isfile(fileList):
    print("Input file fileList " + fileList + " not found")
    sys.exit(-1)


if f1 is not None and f2 is not None:
    m1 = Motif(fileIn=f1)
    m2 = Motif(fileIn=f2)
    fileOut = None
    if out is not None:
        fileOut = open(out, 'w')
    else:
        verbose = True
    m1.fisim(m2, fileOut=fileOut, verbose=verbose)

elif fileList is not None:
    cores = []

    fileList = open(fileList)
    motifs = []
    while not utils.EOF(fileList):
        m = Motif(fileIn=fileList)  # create nas instance of the Motif
        motifs.append(m)
        if core:
            cores.append(m.core(n=5))
    utilsMatrix.createSimMatrix(motifs, cores, out, ID=ID, verbose=verbose)

else:
    helpProgram()
    sys.exit(-1)

