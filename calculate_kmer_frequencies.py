#!/bin/python

#------------------- Description & Notes --------------------#

'''
Description:
    Given a list of FASTA/FASTQ or SAM/BAM files, output (to file) the
    frequencies for each Kmer sequence of length K in a sequence/s.
    Frequencies are calculated either:
        * Across all sequence records (combined).
        * For each sequence record (split).
        * Along the length of a given sequence record (sequential).

Args:
    fastXFiles (filepath):
        List containing the filepath of each file. Files can be
        compressed (.gz) and should contain at least one
        sequence record (FASTA/FASTQ or SAM/BAM). Our processing limit
        seems to be around 3.5 million (3,500,000) sequence records.

    kmerLength (int):
        Length of Kmer sequences. Must be a positive integer. Ideally,
        this should be <= 13 since the total number of Kmer sequences
        exponentially increases (4^K).
            * 4^13 Kmers =    67,108,864  ## Possible
            * 4^14 Kmers =   268,435,456  ## Sometimes possible
            * 4^15 Kmers = 1,073,741,824  ## Probably not possible

Returns:
    oFile (dir):
        Directory containing a list of files. Each file is compressed
        in Parquet format and contains the frequencies for each
        Kmer sequence of length K in a sequence/s.
'''

'''
Design considerations
    Two general approaches:
        * Create table of all possible Kmers, iterate through each Kmer and
          count its frequency in the sequence.
        * Sliding window across the sequence and count the occurrence of
          each Kmer.

    Problems encountered:
        * Table of Kmers does not scale well (4^K Kmers).
        * Sliding window encounters unexpected characters (N, M, etc...)
        * Frequency tables can be big (4^K rows).
'''

'''
Tool evaluation:
    Jellyfish
        * Seems to be the literature standard.
        * Fast, inexact (approximate) Kmer counting.

    KmerCounter (implemented by Brendan's brother)
        * Fast, exact Kmer counting.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import itertools
import sys

# External imports

# Internal imports
from src.kmer import count as kmercount
from src.util import params
from src.util import spark
from src import bio
from src import io

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

def readFiles(fastXFiles, countRevComp):
    seqRecs = io.fastx.read(*fastXFiles)

    if (countRevComp == 'create'):
        seqRecs = (bio.createRevComp(s) for s in seqRecs)
        seqRecs = itertools.chain(*seqRecs)

    elif (countRevComp == 'append'):
        seqRecs = (bio.appendRevComp(s) for s in seqRecs)

    return seqRecs

def combineKmerFrequencies(iFiles, kmerLength, oFile, ignoreNs, countExp):
    print("Combining counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Read FASTX files
            seqRecRdd = sc.parallelize(iFiles)
            seqRecRdd = seqRecRdd.flatMap(io.fastx.read)
            seqRecRdd = kmercount.combined.setup(seqRecRdd)

            ## Get a table containing the kmer counts across all records
            kmerDf    = kmercount.combined.getCounts(seqRecRdd, kmerLength,
                ignoreNs, countExp)
            kmerDf    = kmercount.combined.cleanup(kmerDf, kmerLength)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def splitKmerFrequencies(iFiles, kmerLength, oFile, ignoreNs, countExp):
    print("Splitting counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Read FASTX files
            seqRecRdd = sc.parallelize(iFiles)
            seqRecRdd = seqRecRdd.flatMap(io.fastx.read)
            seqRecRdd = kmercount.split.setup(seqRecRdd, kmerLength)

            ## Get a table containing the kmer counts across for each record
            kmerDf    = kmercount.split.getCounts(seqRecRdd, kmerLength,
                ignoreNs, countExp)
            kmerDf    = kmercount.split.cleanup(kmerDf, kmerLength)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def sequentialKmerFrequencies(iFiles, kmerLength, oFile,
    fastxId, abSize, nWindows, nSteps, ignoreNs, countExp,):

    print("Sequential counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Read FASTX files
            seqRecRdd = sc.parallelize(iFiles)
            seqRecRdd = seqRecRdd.flatMap(io.fastx.read)
            seqRecRdd = kmercount.split.setup(seqRecRdd, kmerLength)

            ## Look for a sequence record if there is a
            ## specific one we want to query
            if (fastxId is not None):
                f = lambda x: x.id == fastxId
                seqRecRdd = seqRecRdd.filter(f)
                if (seqRecRdd.isEmpty()):
                    raise IndexError("Invalid FASTX ID")

            ## Get a table containing the kmer counts across for each record
            seqRecRdd = kmercount.sequential.getSeqRecs(seqRecRdd, abSize,
                nWindows, nSteps)

            kmerDf    = kmercount.split.getCounts(seqRecRdd, kmerLength,
                ignoreNs, countExp)
            kmerDf    = kmercount.split.cleanup(kmerDf, kmerLength)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def main():
    argParser = params.CalculateKmerArgParser()
    argParser.parse()

    if (argParser.cmd == 'combined'):
        argParser.printArgs()

        # seqRecs = readFiles(argParser.iFiles, , argParser.countRevComp)
        combineKmerFrequencies(argParser.iFiles, argParser.kmerLength,
            argParser.oFile, argParser.ignoreNs, argParser.countExp)

    elif (argParser.cmd == 'split'):
        argParser.printArgs()

        # seqRecs = readFiles(argParser.iFiles, argParser.countRevComp)
        splitKmerFrequencies(argParser.iFiles, argParser.kmerLength,
            argParser.oFile, argParser.ignoreNs, argParser.countExp)

    elif (argParser.cmd == 'sequential'):
        argParser.parseSequentialArgs()
        argParser.printArgs()

        # seqRecs = readFiles(argParser.iFiles, argParser.countRevComp)
        sequentialKmerFrequencies(argParser.iFiles, argParser.kmerLength,
            argParser.oFile, argParser.fastxId, argParser.abSize,
            argParser.nWindows, argParser.nSteps,
            argParser.ignoreNs, argParser.countExp)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
