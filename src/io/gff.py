#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import re
from pathlib import Path

# External imports
import pandas as pd

# Internal imports
from .common import isZFile

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*filepaths, asPdf=True):
    fDbs = (_readFile(f, asPdf) for f in filepaths)
    return fDbs

#------------------- Private Classes & Functions ------------#

def _readFile(filepath, asPdf):
    if (isZFile(filepath)):
        stem    = Path(filepath).stem
        gxfType = _getFormat(stem)

    else:
        gxfType = _getFormat(filepath)

    fDb = _getFeatureDB(filepath, asPdf, gxfType)
    return fDb

def _getFormat(filepath):
    if (filepath.endswith('.gtf')):
        return GTF

    elif (filepath.endswith('.gff') \
          or filepath.endswith('.gff3')):
        return GFF

    else:
        raise NotImplementedError("Unknown GXF file")

def _getFeatureDB(filepath, asPdf, gxfType):
    fDb = None
    if (asPdf):
        ## Create DataFrame of the file
        colNames = ['seqname', 'source', 'feature', 'start',
                    'end', 'score', 'strand', 'frame', 'attribute']
        mDf      = pd.read_csv(filepath, sep='\t', comment='#',
            header=None, names=colNames, low_memory=False)

        ## Parse the attributes column
        aDf = mDf['attribute'].to_list()
        if (gxfType == GFF):
            f   = lambda x: dict(p.split('=') for p in x.split(';'))
            aDf = map(f, aDf)
            aDf = pd.DataFrame.from_dict(aDf)

        else:
            f   = lambda x: [re.match('(.*) "(.*)"', y) for y in x.split('; ')]
            g   = lambda x: dict([y.group(1), y.group(2)] for y in x)
            aDf = map(g, map(f, aDf))
            aDf = pd.DataFrame.from_dict(aDf)

        ## Join the tables and adjust column types
        fDb = pd.concat([mDf, aDf], axis=1)
        fDb['start'] = fDb['start'].astype(int)
        fDb['end']   = fDb['end'].astype(int)

    else:
        import gffutils             ## Requires python 3.5; not 3.7
        fDb = gffutils.create_db(filepath,
            ':memory:', force=True, keep_order=True,
            merge_strategy='merge', sort_attribute_values=True)

        # fDb.update(fDb.create_introns())
    return fDb

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
