#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import os

# External imports
import holoviews as hv

# Internal imports
from .common import createDirIfNone, removeFileIfExists

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def write(filepath, plot):
    outputDir = os.path.dirname(filepath)
    createDirIfNone(outputDir)
    removeFileIfExists(filepath)
    hv.save(plot, filepath)

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
