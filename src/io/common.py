#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import gzip
import shutil
from pathlib import Path

# External imports

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Protected Classes & Functions ------------#

def toZFile(filepath):
    with open(filepath, 'rb') as f_input:
        with gzip.open("{}.gz".format(filepath), 'wb') as f_output:
            shutil.copyfileobj(f_input, f_output)

            # Delete the original file after the gzip is done
            Path(filepath).unlink()

def isZFile(filepath):
    if (filepath.endswith('.gz') \
        or filepath.endswith('.gzip')):
        return True

    return False

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
