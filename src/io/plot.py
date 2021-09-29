#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
from pathlib import Path

# External imports
import holoviews as hv

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def write(filepath, plot):
    ## Create DIR & FILE they don't exist
    p = Path(filepath)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.unlink() if p.exists() else None

    hv.save(plot, filepath)

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
