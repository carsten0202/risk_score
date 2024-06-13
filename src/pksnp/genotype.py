###########################################################
#
# ---%%%  allele.py: Module to hold hashable alleles for matching PRS scores  %%%---
#

import logging
import re

from . import Allele

logger = logging.getLogger(__name__)

class Genotype:
	"""
    Class to hold a diploid genotype (ie. two alleles). Can be used as key in e.g. dict(key) = (gt or dosage) or dict(key) = risk_weight
	Must be immutable; do not change after __init__

	genotype (tuple): Tuple of nucleotides, sorted alphabetically. Used for speed in __eq__ and __hash__
    _alleles (tuple): Tuple of alleles
	"""
	def __init__(self, chromosome=None, position=None, rsid=None, genotype=None, build=None):
		"""
	    chromosome (str): The chromosome
	    position (int): Chromosomal position
    	rsid (str): The rsid of the variant
    	genotype (str): A string that can be split into nucleotides; e.g. "A/G"
	    build (str): Genomic build (see _validate)
		"""
		self.genotype = _sort_and_uppercase(re.split("[|/>]", genotype))
		self._alleles = (
			Allele(chromosome=chromosome, position=position, rsid=rsid, allele=self.genotype[0], build=build),
			Allele(chromosome=chromosome, position=position, rsid=rsid, allele=self.genotype[1], build=build),
		)
		Genotype._validate(self)

	def __eq__(self, other):
		if isinstance(other, Genotype):
			return self.alleles == other.alleles
		return NotImplemented

	def __hash__(self):
		return hash(self.alleles)
	
	def __repr__(self):
		return f"Genotype({self.rsid}, {'/'.join(self.genotype)})"

	@staticmethod
	def _validate(self):
		"""What to check?"""
		assert len(self.genotype) == 2, f"ERROR! '{self.genotype}' does not translate to a diploid Genotype."

	@property
	def alleles(self):
		return self._alleles

	@property
	def homozygote(self):
		"""Is the genotype homozygotic (True/False)?"""
		return self._alleles[0] == self._alleles[1]

	@property
	def rsid(self):
		return self._alleles[0].rsid


def _sort_and_uppercase(input_iter):
	"""
	Takes a string or list of strings, converts each to uppercase, sorts them alphabetically, and returns as a tuple.

	Parameters:
	input_iter (iter): String or list of strings to process.

	Returns:
	tuple: Alphabetically sorted tuple of uppercase strings.
	"""
	# Convert each string to uppercase
	uppercased_list = [s.upper() for s in input_iter]
	# Sort the list alphabetically
	sorted_list = sorted(uppercased_list)
	# Convert the list to a tuple and return
	return tuple(sorted_list)
