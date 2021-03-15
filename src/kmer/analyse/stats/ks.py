#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import itertools
import traceback

# External imports
import numpy as np
import pandas as pd
import scipy.stats as sps
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
from .common import *
from ...constants import *
from ... import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPvalues(kmerDf, method='average'):
    from statsmodels.stats import multitest

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Load cluster IDs into memory
    cIds = [x[0] for x in kmerDf.select('ClusterId').distinct().collect()]
    print("NUM CLUSTERS\t{}".format(len(cIds)))

    if (method == 'average'):
        ## Calculate the 'average' Kmer frequencies
        ## This will load the counts into memory
        avgKmerCount = getTotalProbabilities(kmerDf)
        print(avgKmerCount)

        ## Perform Kolmogorov-Smirnov test between each cluster and the 'average'
        print("Running Kolmogorov-Smirnov test between average")
        f     = lambda x, y: _getPvaluesAgainstAverage(x, y, avgKmerCount)
        s     = "ClusterId_x string, ClusterId_y string, p_value double"
        df    = kmerDf.groupby('ClusterId').applyInPandas(f, schema=s)
        hpvDf = df.toPandas()

    elif (method == 'cluster'):
        ## Perform Kolmogorov-Smirnov test between each cluster
        ## and every other cluster.
        print("Running Kolmogorov-Smirnov test between clusters")
        raise NotImplementedError("TODO")

    ## Correct p-values due to multiple testing
    pvalues = hpvDf['p_value'].tolist()
    results = multitest.multipletests(pvalues, alpha=0.01, method='fdr_bh')
    pvalues = results[1]
    hpvDf['p_value'] = pvalues
    hpvDf = hpvDf.sort_values(by=['p_value']).reset_index(drop=True)
    return hpvDf

#------------------- Protected Classes & Functions ------------#

def _getPvaluesAgainstAverage(key, kmerPdf, avgKmerCount):
    print("{}".format(key))

    ## Get the Kmer frequencies from the Pandas DataFrame
    kmerPdf   = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, FILE_COL_NAME,
                         KMER_COL_NAME, COUNT_COL_NAME]]
    kmerPdf   = transform.rotatePdf(kmerPdf)
    kmerPdf   = transform.splitPdf(kmerPdf)
    kmerId    = kmerPdf[0]
    kmerCount = kmerPdf[1]

    ## Calculate the 'average' kmer frequencies of the cluster
    cAvgKmerCount = kmerCount.sum() / kmerId[TOTAL_COL_NAME].sum()
    cAvgKmerCount = cAvgKmerCount.to_frame().T
    cAvgKmerCount.index = [str(key[0])]

    ## Test whether two samples come from a population
    ## with the same distribution (Kolmogorov-Smirnov test)
    ## The null hypothesis is that the two dataset values
    ## are from the distribution
    cDf = pd.concat([cAvgKmerCount, avgKmerCount]).fillna(0)
    kmerCount1 = cDf.loc[cAvgKmerCount.index, :].to_numpy()[0]
    kmerCount2 = cDf.loc[avgKmerCount.index, :].to_numpy()[0]
    r   = sps.ks_2samp(kmerCount1, kmerCount2)
    p   = r[1]

    row = {'ClusterId_x':key, 'ClusterId_y':'Average', 'p_value':p}
    df  = pd.DataFrame(row, index=[0])
    return df

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
