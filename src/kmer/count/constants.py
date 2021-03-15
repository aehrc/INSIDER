#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports

# Internal imports
from ..constants import *

#------------------- Constants ------------------------------#

## Maximum number of output files
MAX_N_FILES = 16

## Maximum number of FASTA/FASTQ entries per Partition
CHUNK_SIZE = 1000000

## Number of Partitions multiplier for Spark
N_PARTS_MULTI = 6

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
