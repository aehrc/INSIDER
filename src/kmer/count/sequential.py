#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np

# Internal imports
from ... import bio

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getSeqRecs(seqRecRdd, abSize, nWindows, nSteps):
    ## Create a (SeqRec) for each (S, E) pair
    ## (SeqRec) => (SeqRec, SeqLen)
    ##          => (SeqRec, (SPos, EPos))
    ##          => (SeqRec)
    f = lambda x: (x, len(x.seq))
    g = lambda x: (getPositions(x, abSize, nWindows, nSteps))
    h = lambda x: bio.fixed.createSpecificSeqRecord(x[0], x[1])
    seqRecRdd = seqRecRdd.map(f).flatMapValues(g).map(h)
    return seqRecRdd

#------------------- Private Classes & Functions ------------#

def getPositions(seqLen, abSize, nWindows, nSteps):
    if (abSize):
        positions  = [(int(round(s)), int(round(s + nWindows)))
                      for s in np.arange(0, seqLen, nSteps)]

    else:
        windowSize = seqLen / nWindows
        stepSize   = windowSize / nSteps
        nSeqPairs  = seqLen - windowSize + 1
        positions  = [(int(round(s)), int(round(s + windowSize)))
                      for s in np.arange(0, nSeqPairs, stepSize)]

    return positions

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
