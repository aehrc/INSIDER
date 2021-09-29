#!/bin/python

#------------------- Description & Notes --------------------#

'''
The Chi-squared test of independence tests the hypothesis that
the observed COUNTS of two (or more) samples are consistent with
the expected COUNTS. Therefore, p-value < 0.05 == reject the hypothesis
that the observed discrepency between samples is due to chance
and hence the samples are different.
'''

'''
Because of the high dimensionality, tiny differences
can lead to statistically significant results (i.e., near-zero
p-values) regardless of the statistical test used. Thus, relying on p-values
for statistical significance may not correlate with the practical
significance anymore.

Instead, we must also quantify the practical significance of the results by
looking at the effect size (e.g., difference between two means).
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import traceback

# External imports
import numpy as np
import scipy.stats as sps
import pyspark.sql.functions as sparkF
import pyspark.sql.types as sparkT
from pyspark.sql import SparkSession

# Internal imports
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPvalues(kmerDfX, kmerDfY=None):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Find all pairs depending on the DataFrames available
    kmerDf = getPairs(kmerDfX, kmerDfY)

    ## Perform Chi-Squared test for each pair
    f = lambda x, y: list(_chi2_n_cramersV(x, y))
    f = sparkF.udf(f, sparkT.ArrayType(sparkT.FloatType()))
    kmerDf = kmerDf.withColumn('r', f('l.count', 'r.count'))
    kmerStatDf = kmerDf.select(sparkF.col('l.seqId').alias('seqId_1'),
        sparkF.col('r.seqId').alias('seqId_2'),
        kmerDf.r[0].alias('p_value'),
        kmerDf.r[1].alias('cramers_v'))
    return kmerStatDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

def _chi2_n_cramersV(x, y):
    try:
        ## Calculate Chi2 for p-value and
        ## Cramer's V (with bias correction) for effect size
        cDf = np.array([x, y])
        (chi2, p) = _chi2(cDf)
        vcorr = _cramersV(cDf, chi2)

    except Exception as e:
        ## Raise any errors that may occur so that we can debug later
        ## As placeholders, we'll set p-value == 1 and vcorr == 0
        # print("Unable to compute, result might be inaccurate")
        # print(e)
        # print(cDf)
        # traceback.print_exc()
        p     = 1
        vcorr = 0

    return (float(p), float(vcorr))

def _chi2(cDf):
    r    = sps.chi2_contingency(cDf)
    chi2 = r[0]
    p    = r[1]
    return (chi2, p)

def _cramersV(cDf, chi2):
    n    = np.sum(cDf)
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

# from statsmodels.stats import multitest
# ## Correct p-values due to multiple testing
# pvalues = hpvDf['p_value'].tolist()
# # results = multitest.multipletests(pvalues, alpha=0.01, method='fdr_bh')
# # pvalues = results[1]
# hpvDf['p_value'] = pvalues
# hpvDf = hpvDf.sort_values(by=['p_value'], ignore_index=True)
# return hpvDf
