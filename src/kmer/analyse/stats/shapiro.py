#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
import scipy.stats as sps
from pyspark.sql import SparkSession

# Internal imports
from .common import *
from ...constants import *
from ... import transform

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getPvalues(kmerDf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## **********
    ## The Shapiro-Wilk test helps us to determine whether one or
    ## more variables are non-normal. Specifically, it tests the hypothesis
    ## that the variables came from a normal distribution. Therefore, p-value
    ## < 0.05 == reject the hypothesis that the variables came from a
    ## normal distribution and hence, the variable is non-normal.
    ##
    ## IMPORTANT: The Shapiro-Wilk test DOES NOT tell us which variables are
    ## non-normal. Other statistical tests are required for us to help
    ## determine which variables are non-normal (i.e., Anderson-Darling,
    ## Pearson chi-square, Kolmogorov-Smirnov)
    ## **********

    ## Perform Shapiro-Wilks test on each group
    ## Although the Shapiro-Wilks and MANOVA tests are mostly sanity checks
    ## we could use the Shapiro-Wilks test could be used as a way to
    ## filter out 'bad' clusters. However, preliminary analyses suggest
    ## that our clusters are always non-normal (but I could be doing
    ## something wrong), so let's be careful with our interpretation
    print("Running Shapiro-Wilks test")
    s     = "ClusterId string, p_value double"
    df    = kmerDf.groupby('ClusterId') \
                  .applyInPandas(_shapiroWilks, schema=s)
    spvDf = df.toPandas()
    return spvDf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

def _shapiroWilks(key, kmerPdf):
    ## Get the Kmer frequencies from the Pandas DataFrame
    kmerPdf   = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, FILE_COL_NAME,
                         KMER_COL_NAME, COUNT_COL_NAME]]
    kmerPdf   = transform.rotatePdf(kmerPdf)
    kmerPdf   = transform.splitPdf(kmerPdf)
    kmerCount = kmerPdf[1]
    kmerCount.index = kmerPdf[0][SEQID_COL_NAME]

    ## Cannot & shouldn't perform test if n <= 3
    if (kmerCount.shape[0] <= 3):
        p = 1

    else:
        ## I think the Shapiro-Wilks in Scipy is only designed for univariate
        ## normality, and not multivariate normality. This may explain
        ## why our p-values are usually 0 and that the results of this
        ## may not be the most informative.
        r = sps.shapiro(kmerCount)
        p = r[1]

    row = {'ClusterId':key, 'p_value':p}
    df  = pd.DataFrame(row, index=[0])
    return df

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
