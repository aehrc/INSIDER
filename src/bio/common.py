#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import random

# External imports
from Bio.Seq import Seq

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def createRevComp(seqRec):
    rcSeqRec  = seqRec.reverse_complement(id=seqRec.id + ':' + 'rc')
    seqRec.id = seqRec.id
    return (seqRec, rcSeqRec)

def appendRevComp(seqRec):
    rcSeqRec     = seqRec.reverse_complement()
    newSeqRec    = seqRec + rcSeqRec
    newSeqRec.id = seqRec.id
    return newSeqRec

def getSeqRec(seqRecs, fastxId=None):
    ## Return a random record if no ID was provided
    if (fastxId is None):
        randFaIdx = random.randint(0, len(seqRecs) - 1)
        seqRec    = [seqRecs[randFaIdx]]
        return seqRec[0]

    else:
        f      = lambda x: x.id == fastxId
        seqRec = list(filter(f, seqRecs))
        if (len(seqRec) == 0):
            return None

        return seqRec[0]

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
