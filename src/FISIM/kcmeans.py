from __future__ import print_function
from __future__ import absolute_import
import sys
import clustering as C


def helpProgram():
    print("Correct usage:\n\tpython kcmeans.py -fileIn <fileSim> -o <fileOut> -k <numCluster>")
    print("\tfileSim is a file containing a similarity matrix with the format of the generated by fisim.py")
    print("\tpython pKcmeans -h prints this help")


if len(sys.argv) > 10 or len(sys.argv) < 3:
    helpProgram()
    sys.exit(-1)

fileIn = None
fileOut = None
verbose = False

i = 1
while i < len(sys.argv):
    if sys.argv[i] == "-fileIn":
        i += 1
        fileIn = sys.argv[i]
    elif sys.argv[i] == "-o":
        i += 1
        fileOut = sys.argv[i]
    elif sys.argv[i] == "-k":
        i += 1
        k = int(sys.argv[i])
    else:
        print("Unknown flag " + sys.argv[i])
        helpProgram()
        sys.exit(-1)
    i += 1

if k < 2:
    print("k should be greater than 2.")
    sys.exit(-1)

C.script_kcmeans(fileIn, k, fileOut)
