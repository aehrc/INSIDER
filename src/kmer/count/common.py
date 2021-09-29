#!/bin/python

#------------------- Description & Notes --------------------#

'''
Kmer frequencies only need to be calcuated in either the forward, or
the reverse direction since they basically contain the same information.
Specifically, for odd K-mers,  Frequencies(forward) == Frequencies(reverse).
* Total N(AAA) = N(AAA, Forward) + N(AAA, Reverse)
               = N(AAA, Forward) + N(TTT, Reverse)

However, for even K-mers, Frequencies(forward) != Frequencies(reverse) because
even K-mers can generate palindromes (Odd K-mers never generate palindromes)
* Total N(AT) = N(AT, Forward) + N(AT, Reverse)
              = N(AT, Forward) + N(AT, Forward)
              = N(AT, Forward) * 2
* Total N(AGCT) = N(AGCT, Forward) + N(AGCT, Reverse)
*               = N(AGCT, Forward) + N(AGCT, Forward)
*               = N(AGCT, Forward) * 2
Thus, for these K-mers, the frequencies are doubled instead of additive.

The above suggests that only odd K-mers should be considered to avoid the
frequency biases due to different strands.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
import pyspark.sql.functions as F
import pyspark.sql.types as T

# Internal imports
from ..common import *

#------------------- Constants ------------------------------#

## Maximum number of output files
MAX_N_FILES = 16

## Maximum number of FASTA/FASTQ entries per Partition
CHUNK_SIZE = 1000000

## Number of Partitions multiplier for Spark
N_PARTS_MULTI = 6

#------------------- Public Classes & Functions -------------#

#------------------- Protected Classes & Functions ------------#

def getObservedSequences(seq, kmerLength):

    """
    Description:
        Generates a list of Kmer sequences of length K in a sequence.

    Args:
        seq (str):
            Sequence to be examined.

        kmerLength (int):
            Length of Kmers. Must be a positive integer.

    Returns:
        kmerSeqs (list):
            List of str containing the Kmer sequences of length K
            in the sequence record
    """

    obsTotal = getObservedTotal(seq, kmerLength)
    return (seq[i:i + kmerLength] for i in range(0, obsTotal))

def getObservedTotal(seq, kmerLength):
    return len(seq) - kmerLength + 1

def getValidRows(kmerSdf):

    """
    Description:
        Generates a Spark DataFrame containing only valid rows. Rows are
        considered invalid if the Kmer sequence contains unexpected characters.

    Args:
        kmerSdf (pyspark.DataFrame):
            Spark DataFrame containing the counts / percentages
            of Kmer sequences.

    Returns:
        kmerSdf (pyspark.DataFrame):
            Filtered Spark DataFrame containing the counts / percentages
            of Kmer sequences.
    """

    ## Ensure that we don't have Kmers containing ambiguous bases
    ## i.e., N's, R's, S's, etc...
    nStr = ''.join(NUCLEOTIDES.keys())
    p    = '^[' + nStr + ']+$'
    return kmerSdf.filter(kmerSdf.kmer.rlike(p))

def getExpectedCounts(kImerDf, kIImerDf, kImerLength, kIImerLength):
    ## Join the tables
    cond   = (F.substring(kImerDf.kImer, 2, kIImerLength) == kIImerDf.kIImer)
    kmerDf = kImerDf.join(kIImerDf, on=[SEQID_COL_NAME, SEQLEN_COL_NAME], how='full') \
        .filter(cond) \
        .withColumnRenamed('kImer', 'preKmer').withColumnRenamed('kImerCount', 'preKmerCount') \
        .withColumnRenamed('kIImer', 'inKmer').withColumnRenamed('kIImerCount', 'inKmerCount')

    cond   = (F.substring(kImerDf.kImer, 1, kIImerLength) == kmerDf.inKmer)
    kmerDf = kmerDf.join(kImerDf, on=[SEQID_COL_NAME, SEQLEN_COL_NAME], how='full') \
        .filter(cond) \
        .withColumnRenamed('kImer', 'suffKmer').withColumnRenamed('kImerCount', 'suffKmerCount')

    ## Calculate the Expected counts
    f = (F.col('preKmerCount') * F.col('suffKmerCount')) / F.col('inKmerCount')
    kmerDf = kmerDf.withColumn(COUNT_COL_NAME, f)
    kmerDf = kmerDf.withColumn(KMER_COL_NAME,
        F.concat(F.substring(kmerDf.preKmer, 1, 1),
            kmerDf.inKmer, F.substring(kmerDf.suffKmer, kImerLength, 1)))

    ## Return the columns we want
    colNames = [SEQID_COL_NAME, SEQLEN_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME]
    kmerDf   = kmerDf.select(colNames)
    return kmerDf

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
