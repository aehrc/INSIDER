#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

REV_COMP = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

#------------------- Public Classes & Functions -------------#

def partitionPdfByRows(pdf, numRows):
    return (pdf.iloc[i:i+numRows, :].reset_index(drop=True)
            for i in range(0, len(pdf), numRows))

def getReverseComplement(seq):
    bases  = reversed(list(seq))
    bases  = [REV_COMP[b] for b in bases]
    revSeq = ''.join(bases)
    return revSeq

def getOuterKmers(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Find Kmers that are in some but not all sequences.
    ## This is effectively an outer join, but we do this by counting
    ## the number of sequences containing each Kmer.
    ## If the count != total number of sequences, then we've found a Kmer
    ## that isn't in every sequence.
    nSeqs  = kmerSdf.select(SEQID_COL_NAME).distinct().count()
    oKmers = kmerSdf.groupby(KMER_COL_NAME).count() \
                    .filter(sparkF.col('count') != nSeqs) \
                    .select(KMER_COL_NAME)
    oKmers = oKmers.rdd.flatMap(lambda x: x).collect()
    return oKmers

def addZeroCounts(key, kmerPdf, oKmers):
    ## Find zero count Kmers using set difference
    currKmers = kmerPdf[KMER_COL_NAME].unique()
    zcKmers   = set(oKmers).difference(currKmers)

    ## Construct dataframe of zero count Kmers
    zcKmers   = pd.DataFrame([zcKmers])
    zcKmers   = zcKmers.T
    zcKmers.columns = [KMER_COL_NAME]
    zcKmers[SEQID_COL_NAME] = str(key[0])
    zcKmers[COUNT_COL_NAME] = 0
    zcKmers[TOTAL_COL_NAME] = kmerPdf[TOTAL_COL_NAME].unique()[0]
    zcKmers[FILE_COL_NAME]  = kmerPdf[FILE_COL_NAME].unique()[0]

    ## Union the tables
    kmerPdf = pd.concat([kmerPdf, zcKmers])
    return kmerPdf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
