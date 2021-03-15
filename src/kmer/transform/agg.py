#!/bin/python

#------------------- Description & Notes --------------------#

'''
We can reduce the number of K-mers by:
* Clustering reverse complement counts (For relatively short K-mers, i.e. K < 5)

We can reduce the number of samples by:
* Clustering counts of related sequences
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys

# External imports
import pandas as pd
import pyspark.sql.functions as sparkF
from pyspark.sql import SparkSession

# Internal imports
from .convert import countsToProbabilities
from ..common import *
from ..constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def countsToDict(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Must be done in Spark because Pandas doesn't allow maptype columns
    ## (ID, Kmer, Percent, *)
    ##     => (ID, Kmer -> Percent, *)
    ##     => (ID, {Kmer -> Percent}, *)
    gCols = kmerSdf.schema.names
    gCols.remove(KMER_COL_NAME)
    gCols.remove(COUNT_COL_NAME)

    schema  = sparkF.struct([KMER_COL_NAME, COUNT_COL_NAME])
    f       = sparkF.map_from_entries(sparkF.collect_list(schema))
    kmerSdf = kmerSdf.groupBy(*gCols).agg(f.alias(KMER_COL_NAME))
    return kmerSdf

def dictToCounts(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## (ID, {Kmer -> Percent}, *) => (ID, Kmer, Percent, *)
    gCols = kmerSdf.schema.names
    gCols.remove(KMER_COL_NAME)

    f        = sparkF.explode(KMER_COL_NAME)
    colNames = f.alias(KMER_COL_NAME, COUNT_COL_NAME)
    kmerSdf  = kmerSdf.select(*gCols, colNames)
    return kmerSdf

def clusterRevCompCounts(kmerPdf):
    ## Create a column containing the reverse complement of each sequence
    kmerPdf['revComp'] = kmerPdf[KMER_COL_NAME].apply(getReverseComplement)

    ## Sum up the counts of complementary K-mers, and double the counts
    ## for the total frequency of each K-mer.
    f = lambda x: '-'.join(sorted(x))
    kmerPdf[KMER_COL_NAME] = kmerPdf[[KMER_COL_NAME, 'revComp']].apply(f, axis=1)
    kmerPdf = kmerPdf.groupby([*KMERDF_COL_NAMES, KMER_COL_NAME]).sum() * 2
    kmerPdf = kmerPdf.reset_index()

    ## Adjust the total counts
    kmerPdf[TOTAL_COL_NAME] = kmerPdf[COUNT_COL_NAME].sum()
    return kmerPdf

def clusterCountsByColumn(kmerSdf, cIdDf):
    cKmerDf = kmerSdf.join(cIdDf, SEQID_COL_NAME, 'left') \
                     .groupby('ClusterId') \
                     .applyInPandas(sumCountsByColumn, schema=kmerSdf.schema)
    return cKmerDf

def sumCountsByColumn(key, kmerPdf):
    ## Sum the counts for each Kmer
    cols       = [KMER_COL_NAME, COUNT_COL_NAME]
    newKmerPdf = kmerPdf[cols].groupby(KMER_COL_NAME).sum()
    newKmerPdf = newKmerPdf.reset_index()

    ## Format the rest of the columns and adjust the total counts
    newKmerPdf[SEQID_COL_NAME] = str(key[0])
    newKmerPdf[FILE_COL_NAME]  = str(key[0])
    newKmerPdf[TOTAL_COL_NAME] = newKmerPdf[COUNT_COL_NAME].sum()
    return newKmerPdf

def averageCountsByColumn(key, kmerPdf):
    ## Sum the counts for each K-mer across sequences
    ## and divide by the number of sequences
    newKmerPdf = sumCountsByColumn(key, kmerPdf)
    numSeqs    = len(kmerPdf[SEQID_COL_NAME].unique())
    newKmerPdf[COUNT_COL_NAME] = newKmerPdf[COUNT_COL_NAME] / numSeqs
    return newKmerPdf

def sumCounts(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Sum the counts for each Kmer
    kmerSdf = kmerSdf.groupby(KMER_COL_NAME) \
                     .agg(sparkF.sum(COUNT_COL_NAME).alias(COUNT_COL_NAME))

    ## Calculate the total number of Kmers
    total   = kmerSdf.select(COUNT_COL_NAME).groupby().sum().collect()[0][0]

    ## Format the rest of the columns and adjust the total counts
    kmerSdf = kmerSdf.withColumn(SEQID_COL_NAME, sparkF.lit('Sum'))
    kmerSdf = kmerSdf.withColumn(FILE_COL_NAME, sparkF.lit('Sum'))
    kmerSdf = kmerSdf.withColumn(TOTAL_COL_NAME, sparkF.lit(total))
    return kmerSdf

def averageCounts(kmerSdf):
    ss = SparkSession.getActiveSession()
    if (ss is None):
        raise EnvironmentError("Must have an active Spark session")

    ## Sum the counts for each K-mer across all sequences and
    ## divide by the number of sequences
    newKmerSdf = sumCounts(kmerSdf)
    numSeqs    = kmerSdf.select(SEQID_COL_NAME).distinct().count()
    f = countsToProbabilities(COUNT_COL_NAME, sparkF.lit(numSeqs))
    newKmerSdf = newKmerSdf.withColumn(COUNT_COL_NAME, f)
    return newKmerSdf

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
