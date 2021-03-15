#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import itertools

# External imports
import numpy as np
import pandas as pd

# Internal imports
from .. import kmer
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getKmerPca(kmerId, kmerCount, **kwargs):
    ## Perform feature selection & reduction
    ## We have to be careful with feature selection (as it changes our plots a lot!)
    # kmerCounts         = ml.feature.select(kmerId, kmerCounts)
    kmerPca = ml.feature.reduce(kmerCount, **kwargs)
    kmerPca = ml.outlier.detect(kmerPca, **kwargs)

    ## Create the analysis table
    kmerPca = pd.concat([kmerId, kmerPca], axis=1)
    return kmerPca

def filterByFilename(pca, *dirPaths):
    filenames = [d.name for d in dirPaths]
    pattern   = '|'.join(filenames)

    cond = pca[kmer.FILE_COL_NAME].str.contains(pattern)
    pca  = pca.drop(pca[cond].index).reset_index(drop=True)
    return pca

def filterByFeatureRange(kmerPca, fRange_1=(0, 1), fRange_2=(0, 1)):
    ## Assume 2 features only (easier to manage)
    cond = (((kmerPca['Component1'] >= fRange_1[0]) & (kmerPca['Component1'] <= fRange_1[1]))
            & ((kmerPca['Component2'] >= fRange_2[0]) & (kmerPca['Component2'] <= fRange_2[1])))
    kmerPca = kmerPca[cond].copy()
    kmerPca = kmerPca.reset_index(drop=True)
    return kmerPca

def filterByOutlier(kmerPca, oLabel=1):
    cond    = (kmerPca[OLABEL_COL_NAME] == oLabel)
    kmerPca = kmerPca.loc[cond].copy()
    kmerPca = kmerPca.reset_index(drop=True)
    return kmerPca

def getParameterGrid(*params):
    pName = [f[0] for f in params]
    pVals = [f[1] for f in params]
    pVals = list(itertools.product(*pVals))

    params = [dict(zip(pName, v)) for v in pVals]
    return params

def getMLColumns(df, mlColname):
    f        = lambda x: mlColname in x
    colNames = list(filter(f, df.columns))
    dataCols = df.loc[:, colNames]
    return dataCols

#------------------- Private Classes & Functions ------------#

def prototypeAnalysis(tKmerPca, tKmerCount):
    isOutlier = tKmerPca[OLABEL_COL_NAME] == -1
    outliers  = tKmerPca.loc[isOutlier]
    oKmerId   = outliers.iloc[:, 0:2]
    oFreq     = tKmerCount.iloc[outliers.index, :]

    oKmerId = oKmerId.reset_index(drop=True)
    oFreq   = oFreq.reset_index(drop=True)

    (oPca, oPcaColNames) = runPca(oFreq)
    oCols                = [oPca]
    oColNames            = [oPcaColNames]

    (oCluster, oClusterColNames) = runClustering(oPca)
    oCols.append(oCluster)
    oColNames.append(oClusterColNames)

    oKmerPca = formatAnalysis(oKmerId, oCols, oColNames)
    return oKmerPca

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()
#------------------------------------------------------------------------------
