"""
Script to extract the negative bed coordinates given a bed file 
distance downstream

Usage:
    python extractnegative.py <Bed-file> <out-file> <downstream-distance>
"""
import os
import sys

def extractnegative(bd,ng,dist):
    '''
    Extract bed file coordinates of the negative sequences
    '''
    with open(bd) as bed:
        with open(ng, "w") as neg:
            for line in bed:
                spl = line.split()
                wr = "%s\t%i\t%i\n" % (spl[0], int(spl[1])+int(dist), int(spl[2])+int(dist))
                neg.write(wr)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    bd = sys.argv[1]
    ng = sys.argv[2]
    dist = sys.argv[3]
    extractnegative(bd, ng, dist)
