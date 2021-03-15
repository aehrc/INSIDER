#!/bin/python

#------------------- Description & Notes --------------------#

'''
Adapted from:
https://github.com/dionresearch/hotelling/blob/master/hotelling/stats.py

I'm not convinced this test is giving me the results I expect / want.
 * We encounter quite a lots of warnings/errors which prevents us from
   calculating a good p-value.
 * P-values seem a bit strange.

'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import itertools
import traceback

# External imports
import numpy as np
import pandas as pd
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
from .common import *
from ...constants import *
from ... import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPvalues(kmerDf, cIdDf=None, method='ca'):
    from hotelling.stats import hotelling_dict
    from statsmodels.stats import multitest

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## **********
    ## The Hotelling T-squared test is the multivariate equivalent of the
    ## univariate T-test. Like the T-test, the Hotelling T-squared tests helps
    ## us to determine whether TWO samples came from the same population.
    ## Specifically, it tests the hypothesis that the samples came from the
    ## same distribution. Therefore, p-value < 0.05 == reject hypothesis that
    ## the samples came from the same distribution, and hence the samples
    ## came from a different distribution.
    ##
    ## There are cases where the p-value for this test COULD be NaN.
    ## In such cases, it is likely that this test is an inappropriate because
    ## samples do not come from a normal distribution. The Shapiro-Wilk test
    ## will help us to confirm this.
    ## **********

    if (cIdDf is None):
        if (method == 'oa'):
            ## To estimate what the true signature might be, sum the counts
            ## of each K-mer across all non-overlapping sequences and divide
            ## the counts by the total number of K-mers.
            ## This will load the counts into memory.
            eKmerDf    = transform.agg.sumCounts(kmerDf)
            eKmerCount = getCountsPdf(eKmerDf)
            eKmerCount = eKmerCount.divide(eKmerCount.sum(axis=1), axis=0)
            print(eKmerCount)

            ## Perform Hotelling T-squared test between each sequence and the 'average'
            print("Running Hotelling T-squared test between obs and average")
            f     = lambda x, y: _getObsEstPvalues(x, y, eKmerCount)
            s     = "ClusterId_x string, ClusterId_y string, p_value double"
            df    = kmerDf.groupby('id').applyInPandas(f, schema=s)
            hpvDf = df.toPandas()

        else:
            raise NotImplementedError('Invalid parameters')

    else:
        ## Load cluster IDs into memory
        cIds = [x[0] for x in cIdDf.select('ClusterId').distinct().collect()]
        print("NUM CLUSTERS\t{}".format(len(cIds)))

        if (method == 'ca'):
            ## To estimate what the true signature might be, sum the counts
            ## of each K-mer across all non-overlapping sequences and divide
            ## the counts by the total number of K-mers.
            ## This will load the counts into memory.
            eKmerDf    = transform.agg.sumCounts(kmerDf)
            eKmerCount = getCountsPdf(eKmerDf)
            eKmerCount = eKmerCount.divide(eKmerCount.sum(axis=1), axis=0)
            print(eKmerCount)

            ## Perform Hotelling T-squared test between each cluster and the 'average'
            print("Running Hotelling T-squared test between average")
            f     = lambda x, y: _getCluEstPvalues(x, y, eKmerCount)
            s     = "ClusterId_x string, ClusterId_y string, p_value double"
            df    = kmerDf.join(sparkF.broadcast(cIdDf), SEQID_COL_NAME, 'left') \
                          .groupby('ClusterId').applyInPandas(f, schema=s)
            hpvDf = df.toPandas()

        elif (method == 'cc'):
            cIdPairs = itertools.combinations(cIds, 2)
            kmerDf   = kmerDf.join(sparkF.broadcast(cIdDf), SEQID_COL_NAME, 'left')

            ## Perform Hotelling T-squared test between each cluster
            ## and every other cluster.
            print("Running Hotelling T-squared test between clusters")
            f        = lambda x: _getCluCluPvalues(kmerDf, x[0], x[1])
            dfs      = map(f, cIdPairs)
            hpvDf    = pd.concat(dfs)

        else:
            raise NotImplementedError('Invalid parameters')

    ## Correct p-values due to multiple testing
    pvalues = hpvDf['p_value'].tolist()
    results = multitest.multipletests(pvalues, alpha=0.01, method='fdr_bh')
    pvalues = results[1]
    hpvDf['p_value'] = pvalues
    hpvDf = hpvDf.sort_values(by=['p_value']).reset_index(drop=True)
    hpvDf['Test'] = 'Hotelling-T2'
    return hpvDf

#------------------- Protected Classes & Functions ------------#

def _getObsEstPvalues(key, kmerPdf, eKmerCount):
    print("{}".format(key))

    ## Get the Kmer frequencies from the Pandas DataFrame
    kmerPdf   = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, FILE_COL_NAME,
                         KMER_COL_NAME, COUNT_COL_NAME]]
    kmerPdf   = transform.rotatePdf(kmerPdf)
    kmerPdf   = transform.splitPdf(kmerPdf)
    kmerCount = kmerPdf[1]
    kmerCount.index = kmerPdf[0][SEQID_COL_NAME]

    ## Calculate the signature of the sequence dividing the counts by
    ## the total number of K-mers. We also artificially increase
    ## the number of observations so that the test can run properly
    kmerCount  = kmerCount.divide(kmerCount.sum(axis=1), axis=0)
    kmerCount.index = [str(key[0])]
    kmerCount = pd.concat([kmerCount] * (eKmerCount.shape[1] + kmerCount.shape[1]))

    ## Calculate P-value
    p   = _hotellingT2(kmerCount, eKmerCount)
    p   = 1 if p == -np.inf else p
    row = {'ClusterId_x':key, 'ClusterId_y':'Average', 'p_value':p}
    df  = pd.DataFrame(row, index=[0])
    return df

def _getCluEstPvalues(key, kmerPdf, eKmerCount):
    print("{}".format(key))

    ## Get the Kmer frequencies from the Pandas DataFrame
    kmerPdf   = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, FILE_COL_NAME,
                         KMER_COL_NAME, COUNT_COL_NAME]]
    kmerPdf   = transform.rotatePdf(kmerPdf)
    kmerPdf   = transform.splitPdf(kmerPdf)
    kmerId    = kmerPdf[0]
    kmerCount = kmerPdf[1]

    ## Calculate the signature of the cluster by summing the
    ## counts of every non-overlapping sequence, and dividing the counts by
    ## the total number of K-mers. We also artificially increase
    ## the number of observations so that the test can run properly
    kmerCount  = kmerCount.sum().to_frame().T
    kmerCount  = kmerCount.divide(kmerCount.sum(axis=1), axis=0)
    kmerCount.index = [str(key[0])]
    kmerCount = pd.concat([kmerCount] * (eKmerCount.shape[1] + kmerCount.shape[1]))

    ## Calculate one sample Hotelling T-squared test
    p   = _hotellingT2(kmerCount, eKmerCount)
    p   = 1 if p == -np.inf else p
    row = {'ClusterId_x':key, 'ClusterId_y':'Average', 'p_value':p}
    df  = pd.DataFrame(row, index=[0])
    return df

def _getCluCluPvalues(kmerDf, cId1, cId2):
    print("{}\t{}".format(cId1, cId2))

    ## Get the Kmer frequencies for each cluster
    ## This will load the counts into memory
    cKmerSdf1   = kmerDf.filter(sparkF.col('ClusterId') == cId1)
    cKmerCount1 = _getCountsPdf(cKmerSdf1)
    cKmerCount1 = cKmerCount1.sum().to_frame().T
    cKmerCount1 = cKmerCount1.divide(cKmerCount1.sum(axis=1), axis=0)
    cKmerCount1.index = [str(key[0])]
    cKmerCount1 = pd.concat([cKmerCount1] * (cKmerCount1.shape[1] * 2))

    cKmerSdf2   = kmerDf.filter(sparkF.col('ClusterId') == cId2)
    cKmerCount2 = _getCountsPdf(cKmerSdf2)
    cKmerCount2 = cKmerCount2.sum().to_frame().T
    cKmerCount2 = cKmerCount2.divide(cKmerCount2.sum(axis=1), axis=0)
    cKmerCount2.index = [str(key[0])]
    cKmerCount2 = pd.concat([cKmerCount2] * (cKmerCount2.shape[1] * 2))

    ## Calculate P-value in both directions and find the maximum of the two
    p1  = _hotellingT2(cKmerCount1, cKmerCount2)
    p2  = _hotellingT2(cKmerCount2, cKmerCount1)
    p   = max(p1, p2)
    p   = 1 if p == -np.inf else p
    row = {'ClusterId_x':cId1, 'ClusterId_y':cId2, 'p_value':p}
    df  = pd.DataFrame(row, index=[0])
    return df

def _hotellingT2(kmerCount1, kmerCount2):
    ## We're unable to compute the statistic if:
    ## * We are comparing the same samples [not likely to happen, see above]
    ## * The sample doesn't have enough observations
    ## For these cases, we set the p-value to -np.inf (see above)
    if (kmerCount1.equals(kmerCount2)):
        p = -np.inf

    elif (kmerCount1.shape[0] == 1):
        p = -np.inf

    else:
        ## If we have relatively few observations, then we may end up removing
        ## all of our variables because they become 'highly corelated'.
        cDf = pd.concat([kmerCount1, kmerCount2]).fillna(0)
        cDf = transform.filter.removeConstantCounts(cDf)
        # cDf = transform.filter.removeCorrelatedCounts(cDf, corr=0.95)

        ## Check whether we have enough observations compared to
        ## the number of variables
        r, c = cDf.shape
        if (r > c):
            p = _ht2(cDf, kmerCount1.index, kmerCount2.index)

        ## One condition of the standard Hotelling T-squared
        ## test is that the number of observations must be greater
        ## than the number of variables (i.e., n1 + n2 - 2 >= p).
        ## However, because of our data, we often break
        ## this condition which results in NaN p-values.
        ##
        ## To account for the above condition, perform the test several
        ## times using randomly selected variables, calculate the p-value
        ## for each test with multiple testing correction and find the
        ## smallest p-value.
        else:
            print("TO REMOVE")
            print(cDf)
            cDfs = (cDf.sample(r - 3, axis=1) for i in range(0, 100))

            f = lambda x: _ht2(x, kmerCount1.index, kmerCount2.index)
            pvalues = map(f, cDfs)
            results = multitest.multipletests(list(pvalues), alpha=0.01,
                method='fdr_bh', returnsorted=True)
            p = results[1][0]

    return p

def _ht2(cDf, idx1, idx2):
    kmerCount1 = cDf.loc[idx1, :].to_numpy()
    kmerCount2 = cDf.loc[idx2, :].to_numpy()

    ## For one-sample tests, transform the data a bit
    ## so that the function can run properly
    if (kmerCount2.shape[0] == 1):
        kmerCount2 = kmerCount2[0]

    try:
        ## Set NaN p-values to 1 to account for multiple correction
        ## testing.
        r = hotelling_dict(kmerCount1, kmerCount2)
        p = 1 if np.isnan(r['p_value']) else r['p_value']

    except np.linalg.LinAlgError as e:
        ## Set p-value to 1 if we encounter this error.
        p = 1

    except Exception as e:
        ## Raise any other type of errors that occur so that we can debug
        ## But set their p-values to 1
        print("An exception occurred")
        print(e)
        traceback.print_exc()
        print(cDf)
        print(kmerCount1)
        print(kmerCount2)
        p = 1

    return p

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
