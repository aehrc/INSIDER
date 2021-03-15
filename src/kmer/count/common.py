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
import itertools

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

NUCLEOTIDES = {'A':0, 'C':1, 'G':2, 'T':3}

#------------------- Public Classes & Functions -------------#

#------------------- Protected Classes & Functions ------------#

def getExpectedSequences(kmerLength):

    """
    Description:
        Generates a list of all possible Kmer sequences of length K.

    Args:
        kmerLength (int):
            Length of Kmers. Must be a positive integer.

    Returns:
        kmerSeqs (generator):
            List of all possible Kmer sequences.

    Raises:
        ValueError:
            If kmerLength is not a positive integer.
    """

    f = itertools.product(NUCLEOTIDES.keys(), repeat=kmerLength)
    return (''.join(c) for c in f)

def getExpectedTotal(kmerLength):
    return (4 ** kmerLength)

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

def getMCMCounts(kImerPdf, kIImerPdf):
    lkImerPdf = kImerPdf.copy()
    rkImerPdf = kImerPdf.copy()

    lkImerPdf[SEQID_COL_NAME] = lkImerPdf[KMER_COL_NAME].str[1:]
    rkImerPdf[SEQID_COL_NAME] = rkImerPdf[KMER_COL_NAME].str[:-1]
    kIImerPdf[SEQID_COL_NAME] = kIImerPdf[KMER_COL_NAME]

    kmerPdf = lkImerPdf.merge(rkImerPdf, on=SEQID_COL_NAME, how='inner')
    kmerPdf = kmerPdf.merge(kIImerPdf, on=SEQID_COL_NAME, how='inner')
    kmerPdf = kmerPdf.drop(columns=[SEQID_COL_NAME, 'total_x', 'total_y', TOTAL_COL_NAME])
    kmerPdf.columns = ['lKmer', 'lProb', 'rKmer', 'rProb', 'mKmer', 'mProb']

    ## The total number of Kmers doesn't make much sense for expected counts
    ## because:
    ##  * Expected counts are not always integer numbers
    ##  * Expected counts / sum(expected counts) != expected probability
    kmerPdf[SEQID_COL_NAME] = kImerPdf[SEQID_COL_NAME].unique()[0]
    kmerPdf[TOTAL_COL_NAME] = 0
    kmerPdf[KMER_COL_NAME]  = kmerPdf['lKmer'] + kmerPdf['rKmer'].str[-1]
    kmerPdf[COUNT_COL_NAME] = (kmerPdf['lProb'] * kmerPdf['rProb']) / kmerPdf['mProb']
    kmerPdf = kmerPdf[[SEQID_COL_NAME, TOTAL_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME]]
    return kmerPdf

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

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
