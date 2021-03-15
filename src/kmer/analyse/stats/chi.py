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

def getPvalues(kmerDf, method='a'):
    from statsmodels.stats import multitest

    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## **********
    ## The Chi-squared test of independence tests the hypothesis that
    ## the observed COUNTS of two (or more) samples are consistent with
    ## the expected COUNTS. Therefore, p-value < 0.05 == reject the hypothesis
    ## that the observed discrepency between samples is due to chance
    ## and hence the samples are different.
    ## **********
    if (method == 'a'):
        ## To estimate what the true signature might be, sum the counts
        ## of each K-mer across all non-overlapping sequences
        ## This will load the counts into memory.
        eKmerDf    = transform.agg.sumCounts(kmerDf)
        eKmerCount = getCountsPdf(eKmerDf)

        ## Perform Chi-Squared test between each sequence and the 'average'
        print("Running Chi-Squared test between obs and average")
        f     = lambda x, y: _getObsEstPvalues(x, y, eKmerCount)
        s     = "id_x string, id_y string, p_value double, cramers_v double"
        df    = kmerDf.groupby('id').applyInPandas(f, schema=s)
        hpvDf = df.toPandas()

    elif (method == 'o'):
        ids     = kmerDf.select(SEQID_COL_NAME).distinct().collect()
        ids     = [x[0] for x in ids]
        idPairs = itertools.combinations(ids, 2)
        print("NUM IDS\t{}".format(len(cIds)))

        ## Perform Chi-Squared test between each sequence and every other sequence
        print("Running Chi-Squared test between obs")
        f     = lambda x: _getObsObsPvalues(kmerDf, x[0], x[1])
        dfs   = map(f, idPairs)
        hpvDf = pd.concat(dfs)

    else:
        raise NotImplementedError('Invalid parameters')

    ## Correct p-values due to multiple testing
    pvalues = hpvDf['p_value'].tolist()
    # results = multitest.multipletests(pvalues, alpha=0.01, method='fdr_bh')
    # pvalues = results[1]
    hpvDf['p_value'] = pvalues
    hpvDf = hpvDf.sort_values(by=['p_value']).reset_index(drop=True)
    return hpvDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

def _getObsEstPvalues(key, kmerPdf, eKmerCount):
    print("{}".format(key))

    ## Get the Kmer frequencies from the Pandas DataFrame
    kmerPdf   = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, FILE_COL_NAME,
                         KMER_COL_NAME, COUNT_COL_NAME]]
    kmerPdf   = transform.rotatePdf(kmerPdf)
    kmerPdf   = transform.splitPdf(kmerPdf)
    kmerCount = kmerPdf[1]
    kmerCount.index = kmerPdf[0][SEQID_COL_NAME]

    ## Compare the observed counts with the expected counts
    eKmerCount = eKmerCount.divide(eKmerCount.sum(axis=1), axis=0)
    eKmerCount = eKmerCount.mul(kmerCount.sum().sum())

    ## Calculate P-value
    p, v = _chi2_n_cramersV(kmerCount, eKmerCount)
    row  = {'id_x':key, 'id_y':'Average',
            'p_value':p, 'cramers_v':v}
    df   = pd.DataFrame(row, index=[0])
    return df

def _getObsObsPvalues(kmerDf, cId1, cId2):
    print("{}\t{}".format(cId1, cId2))

    ## Get the Kmer frequencies for each sequence
    ## This will load the counts into memory
    cKmerSdf1   = kmerDf.filter(kmerDf.id == cId1)
    cKmerCount1 = getCountsPdf(cKmerSdf1)

    cKmerSdf2   = kmerDf.filter(kmerDf.id == cId2)
    cKmerCount2 = getCountsPdf(cKmerSdf2)

    ## Calculate P-value
    p, v = _chi2_n_cramersV(cKmerCount1, cKmerCount2)
    row  = {'id_x':cId1, 'id_y':cId2,
            'p_value':p, 'cramers_v':v}
    df   = pd.DataFrame(row, index=[0])
    return df

def _chi2_n_cramersV(kmerCount1, kmerCount2):
    ## Merge the tables so that we can normalise the columns
    cDf = pd.concat([kmerCount1, kmerCount2]).fillna(0)

    try:
        ## **********
        ## *** Because of the high dimensionality, tiny differences
        ## *** can lead to statistically significant results (i.e., near-zero
        ## *** p-values) regardless of the statistical test used. Therefore,
        ## *** relying on p-values for statistical significance may not
        ## *** correlate with the practical significance anymore.
        ## ***
        ## *** Therefore, to quantify the practical significance of the results,
        ## *** we also need to look at the effect size (e.g., difference between
        ## *** two means).
        ## **********

        ## Calculate Chi2 for p-value and
        ## Cramer's V (with bias correction) for effect size
        (chi2, p)  = _chi2(cDf)
        vcorr = _cramersV(cDf, chi2)

    except Exception as e:
        ## Raise any errors that may occur so that we can debug later
        ## As placeholders, we'll set p-value == 1 and vcorr == 0
        print("An exception occurred")
        print(e)
        traceback.print_exc()
        print(cDf)
        print(kmerCount1)
        print(kmerCount2)
        p     = 1
        vcorr = 0

    return (p, vcorr)

def _chi2(cDf):
    r    = sps.chi2_contingency(cDf.to_numpy())
    chi2 = r[0]
    p    = r[1]
    return (chi2, p)

def _cramersV(cDf, chi2):
    n    = cDf.sum().sum()
    phi2 = chi2 / n
    k, r = cDf.shape

    phi2corr = max(0, phi2 - (((k - 1) * (r - 1)) / (n - 1)))
    kcorr    = k - (((k - 1) ** 2) / (n - 1))
    rcorr    = r - (((r - 1) ** 2) / (n - 1))
    vcorr    = np.sqrt(phi2corr / min((kcorr - 1), (rcorr - 1)))
    return vcorr

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
