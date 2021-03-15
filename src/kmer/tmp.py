#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def countsToCentralisedCounts(oKmerDf, eKmerDf):
    oKmerId    = oKmerDf[0]
    oKmerCount = oKmerDf[1]
    eKmerCount = eKmerDf[1][oKmerCount.columns]

    ## See Reinert et al. (2009)
    ## Xw = Xw - N
    ##  N = (len(X) - k + 1) * P(w)
    eKmerCount = np.multiply(eKmerCount, oKmerId[TOTAL_COL_NAME][:, np.newaxis])
    kmerCount  = np.subtract(oKmerCount, eKmerCount)
    return kmerCount

def countsToStandardisedCounts(oKmerDf, eKmerDf):
    oKmerId    = oKmerDf[0]
    oKmerCount = oKmerDf[1]
    eKmerCount = eKmerDf[1][oKmerCount.columns]

    ## See Reinert et al. (2009) & Ren et al. (2013)
    ## Xw = (Xw - N) / sqrt(N)
    ##  N = (len(X) - k + 1) * P(w)
    ## Not 100% sure whether this is the correct formula
    N = np.multiply(eKmerCount, oKmerId[kmer.TOTAL_COL_NAME][:, np.newaxis])
    kmerCount = np.divide(np.subtract(oKmerCount, N), np.sqrt(N))
    kmerCount = kmerCount.fillna(0)
    return kmerCount

def observedToExpected(kmerCount, **kwargs):
    method = kwargs.pop('method', 'MCM')

    if (method == 'ZOM'):
        oImerCount = kwargs.pop('oImerCount')
        kmerCount  = _observedToExpectedZOM(kmerCount, oImerCount)

    else:
        raise NotImplementedError("Not implemented")

    return kmerCount

#------------------- Protected Classes & Functions ------------#

def _observedToExpectedZOM(oKmerCount, oImerCount):
    colNames    = oKmerCount.columns.tolist()
    kmerLength  = len(colNames[0])
    eKmerCount  = None

    ## We can only calculate the expected ZOM counts
    ## for K-mers > 1
    if (kmerLength > 1):
        eKmerCount = (_getExpectedZOM(c, oImerCount) for c in colNames)
        eKmerCount = pd.concat(eKmerCount, axis=1)

    return eKmerCount

def _getExpectedZOM(c, oImerCount):
    cols = [oImerCount[x] for x in c]

    f = lambda x,y: x * y
    eKmerCount = functools.reduce(f, cols)
    eKmerCount = eKmerCount.rename(c)
    return eKmerCount

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
