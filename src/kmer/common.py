#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports

# Internal imports

#------------------- Constants ------------------------------#

## Column names
SEQID_COL_NAME   = 'seqId'
SEQLEN_COL_NAME  = 'seqLen'
KMER_COL_NAME    = 'kmer'
COUNT_COL_NAME   = 'count'
FILE_COL_NAME    = 'filename'
CID_COL_NAME     = 'ClusterId'

UID_COL_NAME     = 'uId'
SIGID_COL_NAME   = 'sigId'

KMERDF_SEQINFO_COL_NAMES = [SEQID_COL_NAME, SEQLEN_COL_NAME, FILE_COL_NAME]
KMERDF_COL_NAMES         = [*KMERDF_SEQINFO_COL_NAMES, KMER_COL_NAME, COUNT_COL_NAME]

## Base infomation
NUCLEOTIDES = {'A':0, 'C':1, 'G':2, 'T':3}
REV_COMP    = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

#------------------- Public Classes & Functions -------------#

def getReverseComplement(seq):
    bases  = reversed(list(seq))
    bases  = [REV_COMP[b] for b in bases]
    revSeq = ''.join(bases)
    return revSeq

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

def getExpectedTotal(seq_or_kLen):
    ## Check whether we're dealing with a Kmer sequence
    ## or a Kmer length
    if (isinstance(seq_or_kLen, str)):
        seq = seq_or_kLen
        ## Check whether the Kmer sequence is paired with its reverse
        ## complement (i.e., FWD-REV)
        if ('-' in seq):
            kmer = seq.split('-')[0]
            kLen = len(kmer)
            ## Check whether the length of the Kmer sequence is even or odd
            if (kLen % 2 == 0):
                ## Even K-mer lengths can generate palindromes which
                ## changes the total number of possible K-mers.
                ## Based on formulas reported in:
                ## * Apostolou-Karampelis et al. (2019)
                x = (2 * kLen) - 1
                y = (kLen - 1)
                t = (2 ** x) + (2 ** y)
                return t

            else:
                t = (4 ** kLen) / 2
                return int(t)

        else:
            kLen = len(seq)
            return (4 ** kLen)

    else:
        kLen = seq_or_kLen
        return (4 ** kLen)

#------------------- Protected Classes & Functions ------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
