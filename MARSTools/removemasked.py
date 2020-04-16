"""
removemasked.py removes the sequences that have been repeat masked

Takes as input fasta sequence with repeat masks and writes out a fasta
without repeat masked sequences.

Usage:
    python removemasked.py <input-fasta> <output-fasta>

"""
from __future__ import print_function
import sys


def remove_masked(fasta_in, fasta_out):
    """
    Removes fasta sequences that have been repeat-masked
    """
    with open(fasta_in) as fa_in:
        with open(fasta_out, "w") as fa_out:
            for line in fa_in:
                if "N" in line:
                    continue
                # elif len(line.split()) < 3:
                #     continue
                else:
                    fa_out.write(line)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    input_fa = sys.argv[1]
    output_fa = sys.argv[2]
    remove_masked(input_fa, output_fa)