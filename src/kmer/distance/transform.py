#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def matrixToTable(kmerDist):
    cols     = {0:'distance', 'level_0':'id_x', 'level_1':'id_y'}
    kmerDist = kmerDist.stack().reset_index()
    kmerDist = kmerDist.rename(columns=cols)
    return kmerDist

def matrixToSymmetricMatrix(kmerDist, useMax=True):
    ## Get the maximum value for each mirror pair
    triu = np.triu(kmerDist)
    tril = np.tril(kmerDist).T
    m    = np.maximum(triu, tril) if useMax else np.minimum(triu, tril)

    ## Reconstruct the matrix
    mT   = np.copy(m).T
    np.fill_diagonal(mT, 0)
    m    = m + mT

    kmerDist = pd.DataFrame(m, columns=kmerDist.columns, index=kmerDist.index)
    return kmerDist

def tableToMatrix(kmerDist):
    kmerDist = pd.pivot_table(kmerDist, index='id_x', columns='id_y',
        values='distance')
    kmerDist.columns.name = None
    kmerDist.index.name   = None
    return kmerDist

#------------------- Private Classes & Functions ------------#


#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
