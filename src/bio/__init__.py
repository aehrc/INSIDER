#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports

# Internal imports
from .feature import getRandomRecord
from .feature import createSeqRecord
from .fixed import createSpecificSeqRecord
from .fixed import createRandomSeqRecord
from .hybrid import createRandomSeqRecord
from .hybrid import isValid

from .common import createRevComp
from .common import appendRevComp
from .common import getSeqRec

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
