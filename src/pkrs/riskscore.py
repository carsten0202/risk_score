#################################################
#
# --%%%  CLASS: RiskScore  %%%--
#

# Just like the 're' module; you pre-complie a GRS based on a dict of alleles and weights

import logging
import re

logger = logging.getLogger(__name__)

class RiskScore:
	"""
	An algorithm template object for holding the definition of a weight-based risk score.

	weights (dict): Dict with <class Allele>:<float>
	N (float): The scaling denominator
	"""
	pat_genotype = re.compile("([GATC/]+)[:;]+([GATC/]+)")
	pat_variant = re.compile("[GATC]+")
	pat_rsid = re.compile("rs[\d]+")

	def __init__(self, weights, N=None, *args, **kwargs):
		"""
		weights (dict): Dict with weighted alleles
		N (float): Non-zero denominator to be used as scaling factor.
		*args, **kwargs : Any non-matching arguments are happily discarded. A bit dirty, but saves us a lot of trouble in FromPGS classmethod.
		"""
		self.weights = {key: float(value) for key, value in weights.items()}
		self.N = len(self.weights) if N is None else float(N)
		logger.debug(f" Weighted Variants = {self}")
		self._validate(self)
		
	def __repr__(self):
		return str(self.weights)

	@staticmethod
	def _validate(obj):
		"""Validate the Instance"""
		assert obj.N > 0, f"The denominator ('{obj.N}') for the arithmetric mean must be >0. (Set to '1' to disable)"
		assert isinstance(obj.weights, dict), f"The weights must be given as a dict (Allele:weight)."

	def calc(self, sample_data):
		"""
		Calculates a linear polygenic risk score for a single sample.
		
		Parameters:
		sample_data (dict): Sample data with Allele: dosage(float) or gt(tuple).

		Returns:
		dict: A score?
		"""
		prs_score = 0.0
		for allele, genotypes in sample_data.items():
			# For simplicity, assume genotype is a tuple of alleles (e.g., (0, 1) or (1, 1))
			if allele in self.weights:
				prs_score += self.weights.get(allele, 0) * sum(genotypes)
				logger.debug(f"calc: {allele} found. Current Sum = {prs_score}")
		logger.debug(f"calc: Sample={sample_data.sample}, Allelic Sum = {prs_score}")
		return prs_score / self.N

	@property
	def rsids(self):
		"""Return the ids of the weighted alleles"""
		return {allele.rsid:None for allele in self.weights}

	@classmethod
	def FromPGS(cls, pgs, *args, **kwargs):
		"""
		Build a risk score calculator from a PGScatalog riskscore

		cls:             The Class of PGS to construct
		pgs:             Instance of <class 'pgscatalog.core.lib.scorefiles.ScoringFile'>
		*args, **kwargs: Passed on to cls.__init__()
		"""
		from pksnp import Allele

		haplotype = {}
		interaction = {}
		weights = {}

		for variant in pgs.variants:
			build = getattr(pgs, 'genome_build', None)
			allele = str(variant.effect_allele)
			# Should prob subclass or patch the pgs-thingy and access e.g. variant.is_haplotype or variant.is_interaction...
			if cls.pat_genotype.split(allele)[1:3] or allele in getattr(cls, 'complex_alleles', {}): # Currently suboptimal: should be eg: 'if variant.is_interaction or variant.is_diplotype:
				# Load line as interaction...
				cls.register_interaction(variant, interaction, build=build)
			elif allele in getattr(cls, 'haplotype_alleles', {}): # Currently suboptimal: should be eg: 'if variant.is_haplotype:' 
				# Load line as haplotype...
				cls.register_haplotype(variant, haplotype, build=build)
			elif cls.pat_variant.match(allele):
				# Load line as regular variant...
				weights[Allele(chromosome=variant.chr_name, position=variant.chr_position, rsid=variant.rsID, allele=allele, build=build)] = float(variant.effect_weight)
			else:
				logger.error(f" {pgs.pgs_id} contains an effect allele ('{allele}') not supported by class {cls}")
				logger.error(f" You are likely invoking the script with the wrong command, or using the wrong PGS file?")
				import sys
				sys.exit("Exiting due to errors...")

		rs = cls(haplotype=haplotype, interaction=interaction, weights=weights, *args, **kwargs)
		logger.info(f"FromPGS: Read {len(rs.rsids)} weighted variants from {pgs.pgs_id}")
		return rs
