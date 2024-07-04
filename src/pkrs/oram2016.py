#################################################
#
# --%%%  CLASS: Sharp2019  %%%--
#

import logging
logger = logging.getLogger(__name__)

from . import Interaction

class Oram2016(Interaction):
	"""
	Calculate GRS based on Oram et al 2016
	"""
	complex_alleles = {
		"DR3/DR4-DQ8":      "C/T:T/C",
		"DR3/DR3":          "T/T:T/T",
		"DR4- DQ8/DR4-DQ8": "C/C:C/C", # Note the typo in official PGS000021.txt.gz file. PGS is such a well-maintained resource...
		"DR4-DQ8/DR4-DQ8":  "C/C:C/C", # And add in the corrected diplotype in case some pedantic decides to fix the typo...
		"DR4-DQ8/X":        "C/C:T/C",
		"DR3/X":            "C/T:T/T",
	}


