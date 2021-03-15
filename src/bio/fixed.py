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

def createRandomSeqRecord(seqRec, windowSize, limits=None, isFasta=True):
    seq  = ""
    if (windowSize > len(seqRec)):
        print("Adjusting window size to sequence length.")
        windowSize = len(seqRec) - 1

    if (isFasta):
        ## Get a random sequence that doesn't contain unexpected characters
        while (len(seq) == 0 or seq.count("N") != 0):
            f   = _getRandomPosition
            (randStartPos, endPos) = f(seqRec, windowSize, limits)
            seq = seqRec.seq[randStartPos:endPos]

        seqId  = "{}:{}:{}".format(seqRec.id, str(randStartPos), str((endPos - 1)))
        seqRec = SeqRecord(Seq(str(seq)), id=seqId, description='')

    else:
        f       = createRandomSeqRecord
        seqRec  = f(seqRec, windowSize, limits=limits, isFasta=True)
        qValues = [91] * len(seqRec.seq)    ## Random number for a quality
        seqRec.letter_annotations['phred_quality'] = qValues

    return seqRec

def createSpecificSeqRecord(seqRec, pos):
    (startPos, endPos) = pos
    seqId = "{}:{}:{}".format(seqRec.id, str(startPos), str((endPos - 1)))
    seq   = str(seqRec.seq)[startPos:endPos]
    return SeqRecord(Seq(seq), id=seqId, description='')

#------------------- Private Classes & Functions ------------#

def _getRandomPosition(seqRec, windowSize, limits=None):
    if (limits is not None):
        randStartPos = random.randint(limits[0], limits[1])

    else:
        idx          = len(seqRec.seq) - windowSize - 1
        randStartPos = random.randint(0, idx)

    endPos = randStartPos + windowSize
    return (randStartPos, endPos)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
