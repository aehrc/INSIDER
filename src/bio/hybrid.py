#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import random

# External imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def createRandomSeqRecord(sSeqPair, tSeqPair):
    ## Choose a random position to integrate the source
    ## and adjust the length of the new sequence
    idx          = len(tSeqPair[1]) - 1
    randStartPos = random.randint(0, idx)
    endPos       = randStartPos + len(sSeqPair[1])

    ## Randomly insert the source sequence into the target sequence
    seq   = tSeqPair[1][:randStartPos] + sSeqPair[1] + tSeqPair[1][randStartPos:]
    tId   = "{}:{}:{}".format(tSeqPair[0], str(0), str(len(tSeqPair[1]) - 1))
    sId   = "{}:{}:{}".format(sSeqPair[0], randStartPos, endPos)
    seqId = "{}:{}".format(tId, sId)
    return SeqRecord(Seq(seq), seqId, description='')

def isValid(seqPair):
    if (len(seqPair[0].split(':')) == 6):
        return True
    return False

def getSourcePositions(seqPair, windowSize):
    ## Get the position of the source sequence in the hybrid sequence
    ## (tID:StartPos:EndPos:sID:StartPos:EndPos)
    ##     => ID:StartPos:EndPos
    sId      = seqPair[0].split(':')[3:]

    ## To account for reads that can span across both the
    ## source and target, we adjust the start
    startPos = int(sId[1]) - windowSize
    endPos   = int(sId[2])
    return (startPos, endPos)

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
