#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import zlib

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.types as sparkT
import pyspark.sql.functions as sparkF

# Internal imports
from ..constants import *
from ..common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def countsToPercentages(counts: pd.Series, total: pd.Series) -> pd.Series:
    pcts = countsToProbabilities(counts, total)
    pcts = probabilitiesToPercentages(pcts)
    return pcts

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def countsToProbabilities(counts: pd.Series, total: pd.Series) -> pd.Series:
    return ((counts / total))

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def probabilitiesToPercentages(counts: pd.Series) -> pd.Series:
    return ((counts * 100))

@sparkF.pandas_udf(returnType=sparkT.IntegerType())
def percentagesToCounts(pcts: pd.Series, total: pd.Series) -> pd.Series:
    counts = percentagesToProbabilities(pcts)
    counts = probabilitiesToCounts(counts, total)
    return counts

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def percentagesToProbabilities(counts: pd.Series) -> pd.Series:
    return ((counts / 100))

@sparkF.pandas_udf(returnType=sparkT.IntegerType())
def probabilitiesToCounts(counts: pd.Series, total: pd.Series) -> pd.Series:
    return ((counts * total))

@sparkF.pandas_udf(returnType=sparkT.FloatType())
def kmersToComplexity(kmer: pd.Series) -> pd.Series:
    f = lambda x: (len(zlib.compress(x.encode())) - len(x.encode()))
    kmer = kmer.apply(f)
    return kmer

def countsToCentralisedCounts(oKmerDf, eKmerDf):
    oKmerDf = _joinCounts(oKmerDf, eKmerDf)

    ## See Reinert et al. (2009)
    ## Xw = Xw - N
    ##  N = (len(X) - k + 1) * P(w)
    eCount  = oKmerDf['eCount'] * oKmerDf[TOTAL_COL_NAME]
    oKmerDf[COUNT_COL_NAME] = oKmerDf['oCount'] - eCount
    oKmerDf = oKmerDf[[*KMERDF_COL_NAMES, KMER_COL_NAME, COUNT_COL_NAME]]
    return oKmerDf

def countsToStandardisedCounts(oKmerDf, eKmerDf):
    oKmerDf = _joinCounts(oKmerDf, eKmerDf)

    ## See Reinert et al. (2009) & Ren et al. (2013)
    ## Xw = (Xw - N) / sqrt(N)
    ##  N = (len(X) - k + 1) * P(w)
    ## Not 100% sure whether this is the correct formula
    eCount = oKmerDf['eCount'] * oKmerDf[TOTAL_COL_NAME]
    oKmerDf[COUNT_COL_NAME] = np.divide((oKmerDf['oCount'] - eCount), np.sqrt(N))
    oKmerDf[COUNT_COL_NAME] = oKmerDf[COUNT_COL_NAME].fillna(0)
    oKmerDf = oKmerDf[[*KMERDF_COL_NAMES, KMER_COL_NAME, COUNT_COL_NAME]]
    return kmerCount



#------------------- Protected Classes & Functions ------------#

def observedToExpected(kmerCount, **kwargs):
    method = kwargs.pop('method', 'MCM')

    if (method == 'ZOM'):
        oImerCount = kwargs.pop('oImerCount')
        kmerCount  = _observedToExpectedZOM(kmerCount, oImerCount)

    else:
        raise NotImplementedError("Not implemented")

    return kmerCount

#------------------- Private Classes & Functions ------------#

def _joinCounts(oKmerDf, eKmerDf):
    oKmerDf = oKmerDf.rename(columns={COUNT_COL_NAME:'oCount'})
    eKmerDf = eKmerDf.rename(columns={COUNT_COL_NAME:'eCount'}) \
                     .drop(columns=[TOTAL_COL_NAME, FILE_COL_NAME])

    ## Join tables
    cond    = [SEQID_COL_NAME, KMER_COL_NAME]
    oKmerDf = oKmerDf.merge(eKmerDf, on=cond, how='outer')

    ## Fix up some of the values
    colNames = [TOTAL_COL_NAME, FILE_COL_NAME]
    oKmerDf[colNames] = oKmerDf[colNames].fillna(method='ffill')
    oKmerDf['oCount'] = oKmerDf['oCount'].fillna(0)
    return oKmerDf

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

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
