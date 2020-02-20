from __future__ import print_function
from numpy import *


def EOF(f):
    pos = f.tell()
    print(pos, end=' ')
    if f.readline():
        f.seek(pos)
        return False
    else:
        return True


def decideFSMorNot(mat):
    return sum(mat[0, :]) < 1.5


