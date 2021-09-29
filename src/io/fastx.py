#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import gzip
import itertools
from pathlib import Path

# External imports
from Bio import SeqIO

# Internal imports
from .common import isZFile, toZFile

#------------------- Constants ------------------------------#

FASTA = 'fasta'
FASTQ = 'fastq'

#------------------- Public Classes & Functions -------------#

def read(*filepaths, **kwargs):
    seqRecs = (_readFile(f, **kwargs) for f in filepaths)
    seqRecs = itertools.chain(*seqRecs)
    return seqRecs

def write(filepath, seqRecords, fastXtype, compress=False):
    ## Create DIR & FILE they don't exist
    p = Path(filepath)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.unlink() if p.exists() else None

    ## Write file
    SeqIO.write(seqRecords, filepath, fastXtype)
    if (compress):
        toZFile(filepath)

#------------------- Private Classes & Functions ------------#

def _readFile(filepath, **kwargs):

    """
    Description:
        Reads a FASTX file/s and generates a list of sequence records.

    Args:
        filepath (str):
            Filepath string.

    Returns:
        seqRecs (list<Bio.SeqRecord.SeqRecord>):
            List of sequence records.
    """

    if (isZFile(filepath)):
        stem      = Path(filepath).stem
        fastxType = _getFormat(stem)
        with gzip.open(filepath, 'rt') as filehandle:
            seqIter = SeqIO.parse(filehandle, fastxType, **kwargs)
            yield from seqIter

    else:
        fastxType = _getFormat(filepath)
        seqIter   = SeqIO.parse(filepath, fastxType, **kwargs)
        yield from seqIter

def _getFormat(filepath):
    if (filepath.endswith('.fa') \
        or filepath.endswith('.fasta')):
        return FASTA

    elif (filepath.endswith('.fq') \
          or filepath.endswith('.fastq')):
        return FASTQ

    else:
        raise NotImplementedError("Unknown fastX file")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
