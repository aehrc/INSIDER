#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import gzip
import itertools
from pathlib import Path

# External imports
import pysam

# Internal imports
from .common import isZFile

#------------------- Constants ------------------------------#

BAM = 'b'
SAM = 's'

#------------------- Public Classes & Functions -------------#

def read(*filepaths, **kwargs):
    alignRecs = [_readFile(f, **kwargs) for f in filepaths]
    alignRecs = itertools.chain(*alignRecs)
    return alignRecs

#------------------- Private Classes & Functions ------------#

def _readFile(filepath, **kwargs):
    if (isZFile(filepath)):
        stem    = Path(filepath).stem
        xamType = _getFormat(stem)

    else:
        xamType = _getFormat(filepath)

    alignF    = pysam.AlignmentFile(filepath, 'r' + xamType)
    alignIter = alignF.fetch(until_eof=True, **kwargs)
    yield from alignIter

def _getFormat(filepath):
    if (filepath.endswith('.bam')):
        return BAM

    elif (filepath.endswith('.sam')):
        return ''

    else:
        raise NotImplementedError("Unknown Xam file")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
