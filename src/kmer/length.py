#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import math

# External imports
import pandas as pd
import holoviews as hv
from holoviews import opts
from scipy.stats import entropy

# Internal imports
from .constants import *
from .transform.convert import percentagesToCounts
from .transform.convert import observedToExpected

#------------------- Constants ------------------------------#

LENGTH_COL_NAME      = 'length'
UNIQUE_COL_NAME      = 'numUnique'
NONUNIQUE_COL_NAME   = 'numNonunique'
TOTAL_KMERS_COL_NAME = 'numKmers'
CRE_COL_NAME         = 'cre'
STAT_COL_NAMES       = [UNIQUE_COL_NAME, NONUNIQUE_COL_NAME, TOTAL_KMERS_COL_NAME]

#------------------- Public Classes & Functions -------------#

def selectFromSequenceLength(*seqs):

    """
    Description:
        Select the optimal length of K based on the average sequence length.
        Based on the formulas reported by Luczak et al. (2018).
        and Sims et al. (2009)

    Args:
        seqs (list):
            List of sequences.

    Returns:
        k (int):
            Optimal length of K for a given set of sequences
    """

    seqLens    = [len(s) for s in seqs]
    meanSeqLen = sum(seqLens) / len(seqs)
    k          = math.ceil(math.log(meanSeqLen, 4)) - 1
    return k

def selectRangeFromCounts(kmerIds, kmerCounts, kmerLengths):

    """
    Description:
        Select the length of K range based on the number of non-unique
        Kmers and the cumulative relative entropy (CRE). Based on the formulas
        reported by Sims et al. (2009).

    Args:
        kmerIds (list):
            List of kmerIds.

        kmerCounts (list):
            List of kmerCounts.

        kmerLengths (list):
            List of kmer lengths.

    Returns:
        (ull, lll) (tupple):
            Length of K range for a given set of sequences.
    """

    ull      = getUpperLimit(kmerIds, kmerCounts, kmerLengths)
    lll      = getLowerLengthLimit(kmerIds, kmerCounts, kmerLengths)
    return (ull, lll)

def getUpperLimit(kmerIds, kmerCounts, kmerLengths):
    kmerDict = _toKmerDict(kmerIds, kmerCounts, kmerLengths)
    cre      = _getCRE(kmerDict)

    ## Find the point where CRE == 0, this corresponds to the
    ## upper length limit according to Sims et al. (2009)
    cond = cre['cre'] == 0
    ull  = cre[cond]['length'].unique()
    ull  = max(ull) if len(ull) != 0 else 0
    return ull

def getLowerLimit(kmerIds, kmerCounts, kmerLengths):
    kmerDict = _toKmerDict(kmerIds, kmerCounts, kmerLengths)
    mDfs     = [_getMetrics(kmerDf, k) for k, kmerDf in kmerDict.items()]
    mDf      = pd.concat(mDfs)

    ## Find the point where nNonuniq is maximum, this corresponds
    ## to the lower length limit according to Sims et al. (2009)
    idxs = mDf.groupby(SEQID_COL_NAME)[NONUNIQUE_COL_NAME] \
              .transform(max) == mDf[NONUNIQUE_COL_NAME]
    lll  = mDf[idxs][LENGTH_COL_NAME].mean()
    lll  = math.ceil(lll)
    return lll

def getUpperLimitPlot(kmerIds, kmerCounts, kmerLengths):
    kmerDict = _toKmerDict(kmerIds, kmerCounts, kmerLengths)
    cre      = _getCRE(kmerDict)

    d  = hv.Dataset(cre, LENGTH_COL_NAME, CRE_COL_NAME)
    c  = d.to(hv.Curve, LENGTH_COL_NAME, CRE_COL_NAME, groupby=SEQID_COL_NAME)
    c  = c.overlay()
    return c

def getLowerLimitPlot(kmerIds, kmerCounts, kmerLengths):
    kmerDict = _toKmerDict(kmerIds, kmerCounts, kmerLengths)
    mDfs     = [_getMetrics(kmerDf, k) for k, kmerDf in kmerDict.items()]
    mDf      = pd.concat(mDfs)

    plots = {}
    for m in STAT_COL_NAMES:
        d  = hv.Dataset(mDf, LENGTH_COL_NAME, m)
        c  = d.to(hv.Curve, LENGTH_COL_NAME, m, groupby=SEQID_COL_NAME)
        plots[m] = c.overlay()

    plots = hv.HoloMap(plots, kdims='metrics')
    return plots

#------------------- Private Classes & Functions ------------#

def _toKmerDict(kmerIds, kmerCounts, kmerLengths):
    kmerDf   = [(kmerId, kmerCount)
                for kmerId, kmerCount in zip(kmerIds, kmerCounts)]
    kmerDict = {k:kmerDf for kmerDf, k in zip(kmerDf, kmerLengths)}
    return kmerDict

def _getMetrics(kmerDf, k):
    kmerCount  = percentagesToCounts(kmerDf)
    nUniq      = kmerCount.apply(_countUnique, axis=1)
    nNonuniq   = kmerCount.apply(_countNonunique, axis=1)
    totalKmers = nUniq + nNonuniq

    cCols = pd.concat([nUniq, nNonuniq, totalKmers], axis=1)
    cCols.columns = STAT_COL_NAMES

    idCol = kmerDf[0][[SEQID_COL_NAME]]
    idCol[LENGTH_COL_NAME] = k
    mDf   = pd.concat([idCol, cCols], axis=1)
    return mDf

def _countUnique(kmerRow):
    kmers       = kmerRow.to_list()
    uniqueKmers = filter(lambda x: x == 1, kmers)
    nKmers      = len(list(uniqueKmers))
    return nKmers

def _countNonunique(kmerRow):
    kmers          = kmerRow.to_list()
    nonuniqueKmers = filter(lambda x: x > 1, kmers)
    nKmers         = len(list(nonuniqueKmers))
    return nKmers

def _getCRE(kmerDict):
    re  = [_getRelativeEntropy(kmerDict, k) for k in kmerDict if k > 2]
    re  = pd.concat(re, axis=1)

    cre = re.cumsum(axis=1)
    cre = pd.concat([kmerDict[1][0], re], axis=1)
    cre = pd.melt(cre, id_vars=SEQID_COL_NAME, value_vars=re.columns.tolist(),
                       var_name=LENGTH_COL_NAME, value_name=CRE_COL_NAME)
    return cre

def _getRelativeEntropy(kmerDict, k):
    oKmerCount   = kmerDict[k][1]
    oKImerCount  = kmerDict[k-1][1]
    oKIImerCount = kmerDict[k-2][1]
    eKmerCount   = observedToExpected(oKmerCount, method='MCM',
        oKImerCount=oKImerCount, oKIImerCount=oKIImerCount)
    e = entropy(eKmerCount, oKmerCount, axis=1)
    e = pd.DataFrame(e, columns=[k]) \
          .replace([np.inf, np.nan], 0)
    return e

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
