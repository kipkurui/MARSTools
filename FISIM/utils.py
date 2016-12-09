from numpy import *


def EOF(f):
    pos = f.tell()
    print pos,
    if f.readline():
        f.seek(pos)
        return False
    else:
        return True


def decideFSMorNot(mat):
    return sum(mat[0, :]) < 1.5


