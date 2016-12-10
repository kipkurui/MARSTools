import types

from numpy import *

import fuzzyIntegral

import utils


class Motif(object):

    def __init__(self, fileIn=None, matrix=None, name=None, ID=None):

        self.matFSM = None
        self.matPSSM = None
        self.FSM = None
        self.n = None
        self.numSamples = None
        self.name = name
        self.ID = ID

        if fileIn is not None:
            self.readFromFile(fileIn)
        elif matrix is not None:
            self.FSM = utils.decideFSMorNot(matrix)

            if self.FSM:
                numSamples = sum(matrix[0, :])
                self.matFSM = array(matrix, float)
                self.n = matrix.shape[0]
                self.calculatePSSMFromFSM(numSamples)

            else:
                self.numSamples = sum(matrix[0, :])
                self.matPSSM = array(matrix, int)

                self.n = matrix.shape[0]
                self.calculateFSMFromPSSM()

        else:
            print "Wrong usage for Motif.__init__()"

    def readFromFile(self, fileIn):
        if type(fileIn) == types.StringType:
            fileIn = open(fileIn)
        flag = 0

        while flag == 0:
            line = fileIn.readline()
            if line == "\n":
                line = fileIn.readline()
            #if line.split(" ")[0] == 'URL':
                #line = fileIn.readline()
            if line.split(" ")[0] == "MOTIF":
                header = line.split()
                name = header[1]
                self.ID = name
                flag = 1
                line = fileIn.readline()

        mat = []
        check = 0
        col = line.split()
        while flag == 1:
            line = fileIn.readline()
            col = line.split()
            if line == "\n":
                line = fileIn.readline()
                col = line.split()
            if col[0] == "letter-probability":
                w = col[5]
                flag += 1
                line = fileIn.readline()
                col = line.split()
        while col and (col[0] != 'MOTIF' or col[0] != 'URL' or col[0] != 'MEME'):
            mat.append([float(s) for s in col])
            pos = fileIn.tell()
            if line == "\n":
                fileIn.readline()
            col = fileIn.readline().split()
        mat = array(mat, float)

        self.FSM = utils.decideFSMorNot(mat)

        if self.FSM:
            self.numSamples = sum(mat[0, :])
            self.matFSM = mat.astype(float)
            self.n = mat.shape[0]
            self.calculatePSSMFromFSM()
            self.name = name
        else:
            self.numSamples = sum(mat[0,:])
            self.matPSSM = array(mat, float)
            self.n = mat.shape[0]
            self.calculateFSMFromPSSM()
            self.name = name

    def calculateFSMFromPSSM(self):
        self.matFSM = array(self.matPSSM,float)
        self.matFSM /= self.numSamples

    def calculatePSSMFromFSM(self, numSamples=None):
        if numSamples == None:
            numSamples = self.numSamples

        self.numSamples = numSamples

        self.matPSSM=self.matFSM

    def fisim(self, otherMotif, bigger=False, fileOut = None, verbose = False):
        """Computes the Fuzzy Integral Similarity as defined in Garcia et al (2009)
Returns a list containing:
    list[0] (float) -> Value of FISim.
    list[1] (int) -> Position where the optimal alignment starts.
    list[2] (boolean) -> True if the optimal value was found in the reverse opposite sequence."""
        if (self.matFSM.shape[0] <= otherMotif.matFSM.shape[0]):
            minSize = self.matFSM.shape[0]
            times = otherMotif.matFSM.shape[0] - minSize + 1
        else:
            return otherMotif.fisim(self, True, fileOut, verbose)

        bestRev = False
        curTime = 0
        startPosition = 0
        maxSim = -1

        while (curTime < times):
            curSim = 0.0

            for i in range(minSize):
                importances = [max(value) for value in zip(self.matFSM[i],otherMotif.matFSM[i+curTime])]
                distances = [1-abs(value[0]-value[1]) for value in zip(self.matFSM[i],otherMotif.matFSM[i+curTime])]

                curSim += fuzzyIntegral.fuzzyIntegral(importances, distances)

            if curSim > maxSim:
                maxSim = curSim
                startPosition = curTime

            curTime += 1

        revMatrix = self.calculateRevMatrixFSM()
        curTime = 0

        while curTime < times:
            curSim = 0.0
            for i in range(minSize):

                importances = [max(value) for value in zip(revMatrix[i],otherMotif.matFSM[i+curTime])]
                distances = [1-abs(value[0]-value[1]) for value in zip(revMatrix[i],otherMotif.matFSM[i+curTime])]

                curSim += fuzzyIntegral.fuzzyIntegral(importances, distances)

            if curSim > maxSim:
                bestRev = True
                maxSim = curSim
                startPosition = curTime

            curTime += 1

        if fileOut is not None:
            fileOut.write("FISim value: " + str(maxSim/minSize) + "\n")
            fileOut.write("Start position: " + str(startPosition) + "\n")
            fileOut.write("Reverse opposite: " + str(bestRev) + "\n")

        if verbose:
            print "FISim value: " + str(maxSim/minSize)
            print "Start position: " + str(startPosition)
            print "Reverse opposite: " + str(bestRev)

        return maxSim/minSize, startPosition, bestRev


    def calculateRevMatrixFSM(self):
        ret = array(zeros(4*self.matFSM.shape[0],float))
        ret.shape = (self.matFSM.shape[0],4)
        for i in range(self.matFSM.shape[0]):
                for j in range(4):
                    ret[i, j] = self.matFSM[self.matFSM.shape[0]-1-i,3-j]

        return ret

    def calculateRevMatrixPSSM(self):
        ret = array(zeros(4*self.matPSSM.shape[0],float))
        ret.shape = (self.matPSSM.shape[0],4)
        for i in range(self.matPSSM.shape[0]):
                for j in range(4):
                    ret[i,j] = self.matPSSM[self.matPSSM.shape[0]-1-i,3-j]

        return ret

    def core(self, n=5):
        if n > self.n:
            n = self.n

        coreStart = 0
        bestConservation = 0
        for i in range(self.n - n):
            curConservation = 0.0
            for j in range(n):
                curConservation += max(self.matPSSM[i+j,:])

            if curConservation > bestConservation:
                bestConservation = curConservation
                coreStart = i
        cor = Motif(matrix=self.matPSSM[coreStart:coreStart+n, :], name=self.name+'_core', ID=self.ID+'_core')
        return cor
