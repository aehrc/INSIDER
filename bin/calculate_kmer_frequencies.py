#!/bin/python

#------------------- Description & Notes --------------------#

'''
Description:
    Given a list of FASTA/FASTQ or SAM/BAM files, output (to file) the
    frequencies for each oligonucleotide sequence of length K in a sequence/s.
    Frequencies are calculated either:
        * Across all sequence records (combined).
        * For each sequence record (split).

Args:
    fastXFiles (filepath):
        List containing the filepath of each file. Files can be
        compressed (.gz) and should contain at least one
        sequence record (FASTA/FASTQ or SAM/BAM). Our processing limit
        seems to be around 3.5 million (3,500,000) sequence records.

    kmerLength (int):
        Length of oligonucleotide sequences. Must be a positive integer.
        Ideally, this should be <= 13 since the total number of possible
        oligonucleotide sequences exponentially increases (4^K).
            * 4^13 Kmers =    67,108,864  ## Possible
            * 4^14 Kmers =   268,435,456  ## Sometimes possible
            * 4^15 Kmers = 1,073,741,824  ## Probably not possible

Returns:
    oFile (dir):
        Directory containing a list of files. Each file is compressed
        in Parquet format and contains the frequencies for each
        oligonucleotide sequence of length K in a sequence/s.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import argparse
import sys
from pathlib import Path

# External imports

# Internal imports
THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str(THIS_DIR.parent))       ## Allow us to import from SRC
from src.kmer import count as kmercount
from src.util import spark
from src import io

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def combineKmerFrequencies(iFiles, kmerLength, oFile, ignoreNs, countExp):
    print("Combining counts")
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
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
    params = get_spark_params()
    with spark.getSparkSession(params) as ss:
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

#------------------- Private Classes & Functions ------------#

def get_spark_params():
    # import os                                  ## For HPC
    params = [
        ## Driver
        # ('spark.driver.cores', '5'),           ## Same as executor
        ('spark.driver.memory', '27G'),        ## Same as executor
        ('spark.driver.maxResultSize', '0'),

        ## Executor
        # ('spark.executor.cores', '5'),         ## Didn't seem to have any effect...
        # ('spark.executor.instances', '99'),    ## Didn't seem to have any effect...
        ('spark.executor.memory', '27G'),
        ('spark.executor.heartbeatInterval', '60s'),

        ## SQL
        ('spark.sql.broadcastTimeout', '600s'),
        ('spark.sql.execution.arrow.pyspark.enabled', 'true'),
        ('spark.sql.shuffle.partitions', '200'),

        ## Misc
        ('spark.local.dir', THIS_DIR),
        # ('spark.local.dir', os.environ['MEMDIR']), ## For HPC
        ('spark.network.timeout', '600s'),
        # ('spark.default.parallelism', '8')
    ]
    return params

def make_parser():
    def _initArgs(p):
        p.add_argument("-f", help="FASTA/FASTQ files (.fa/.fq)",
            nargs='+', type=str, required=True)
        p.add_argument("-k", help="Kmer length",
            type=int, required=True)
        p.add_argument("-o", help="Output file (.snappy.parquet)",
            type=str, required=True)
        p.add_argument("-n", help="Ignore Kmers containing ambiguous bases \
            (i.e., N's)", action='store_true')
        p.add_argument("-e", help="Calculate expected Kmer counts instead \
            of observed Kmer counts. Expected counts are based on the (k-2) \
            Markov Chain Model", action='store_true')

    parser    = argparse.ArgumentParser(description='Compute \
        oligonucleotide frequencies of sequence/s')
    subparser = parser.add_subparsers(dest='command')
    combined  = subparser.add_parser('combined')
    split     = subparser.add_parser('split')
    _initArgs(combined)
    _initArgs(split)
    return parser

def main(parser):
    args = parser.parse_args()
    if (args.command is None):
        parser.print_help()
        sys.exit(1)

    else:
        if (args.k < 0):
            parser.print_help()
            parser.error('Invalid value. K > 0.')
            sys.exit(1)

        if (args.e and args.k < 3):
            parser.print_help()
            parser.error('Expected counts can only be calculated for K > 3.')
            sys.exit(1)

        if (args.command == 'combined'):
            print(args)
            combineKmerFrequencies(args.f, args.k, args.o, args.n, args.e)

        elif (args.command == 'split'):
            print(args)
            splitKmerFrequencies(args.f, args.k, args.o, args.n, args.e)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    parser = make_parser()
    main(parser)

#------------------------------------------------------------------------------
