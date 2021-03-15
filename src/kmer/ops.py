#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd

# Internal imports
from .constants import *
from .common import partitionPdfByRows

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getRowsByValue(kmerDf, colName, *values):
    if (len(values) != 0):
        pattern = '|'.join(values)

        ## If kmerDf is a single dataframe
        if (isinstance(kmerDf, pd.DataFrame)):
            cond   = kmerDf[colName].str.contains(pattern)
            kmerDf = kmerDf[cond].reset_index(drop=True)

        ## Otherwise, we're dealing with a
        ## tuple consisting of (KmerId, KmerCount)
        else:
            kmerId    = kmerDf[0]
            kmerCount = kmerDf[1]

            cond      = kmerId[colName].str.contains(pattern).fillna(False)
            kmerId    = kmerId[cond].reset_index(drop=True)
            kmerCount = kmerCount[cond].reset_index(drop=True)
            kmerDf    = (kmerId, kmerCount)

    return kmerDf

def sortRowsByColumn(kmerDf, colName):
    ## If kmerDf is a single dataframe
    if (isinstance(kmerDf, pd.DataFrame)):
        kmerDf = kmerDf.sort_values(by=[colName])
        kmerDf = kmerDf.reset_index(drop=True)

    ## Otherwise, we're dealing with a
    ## tuple consisting of (KmerId, KmerCount)
    else:
        kmerId    = kmerDf[0]
        kmerCount = kmerDf[1]

        kmerId    = kmerId.sort_values(by=[colName])
        kmerCount = kmerCount.reindex(kmerId.index)

        kmerId    = kmerId.reset_index(drop=True)
        kmerCount = kmerCount.reset_index(drop=True)
        kmerDf    = (kmerId, kmerCount)

    return kmerDf

def partitionByRows(kmerDf, numRows=200):
    ## If kmerDf is a single dataframe
    if (isinstance(kmerDf, pd.DataFrame)):
        kmerDfs = partitionPdfByRows(kmerDf, numRows)

    ## Otherwise, we're dealing with a
    ## tuple consisting of (KmerId, KmerCount)
    else:
        kmerId    = kmerDf[0]
        kmerCount = kmerDf[1]

        kmerCounts = partitionPdfByRows(kmerCount, numRows)
        kmerIds    = partitionPdfByRows(kmerId, numRows)
        kmerDfs    = ((kmerId, kmerCount)
                      for kmerId, kmerCount in zip(kmerIds, kmerCounts))

    return kmerDfs

def addRows(*kmerDfs):
    t = type(kmerDfs[0])
    if (not all(isinstance(kmerDf, t) for kmerDf in kmerDfs)):
        raise TypeError('Unable to concat. Not implemented')

    else:
        if (isinstance(kmerDfs[0], pd.DataFrame)):
            kmerdf = pd.concat(kmerDfs) \
                       .drop_duplicates(ignore_index=True)

        else:
            kmerIds    = [kmerDf[0] for kmerDf in kmerDfs]
            kmerCounts = [kmerDf[1] for kmerDf in kmerDfs]

            kmerId     = pd.concat(kmerIds).reset_index(drop=True)
            kmerCount  = pd.concat(kmerCounts).reset_index(drop=True)
            kmerDf     = pd.concat([kmerId, kmerCount], axis=1) \
                           .drop_duplicates(ignore_index=True)

            kmerId     = kmerDf[kmerId.columns]
            kmerCount  = kmerDf[kmerCount.columns]
            kmerDf     = (kmerId, kmerCount)

    return kmerDf

def sampleRows(kmerDf, n=100, random_state=42):
    ## If kmerDf is a single dataframe
    if (isinstance(kmerDf, pd.DataFrame)):
        kmerDf = kmerDf.sample(n=n, random_state=random_state)

    ## Otherwise, we're dealing with a
    ## tuple consisting of (KmerId, KmerCount)
    else:
        kmerId    = kmerDf[0]
        kmerCount = kmerDf[1]

        kmerId    = kmerId.sample(n=n, random_state=random_state)
        kmerCount = kmerCount.sample(n=n, random_state=random_state)

        kmerId    = kmerId.reset_index(drop=True)
        kmerCount = kmerCount.reset_index(drop=True)
        kmerDf    = (kmerId, kmerCount)

    return kmerDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
