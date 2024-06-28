#################################################
#
# --%%%  CLASS: Sharp2019  %%%--
#

import logging
logger = logging.getLogger(__name__)

from . import Interaction

class Sharp2019(Interaction):
	"""
	Calculate GRS based on Sharp et al 2019
	
	haplotype : Haplotype data from Table S1 in Sharp2019 et al. To be used *only* when no interaction is present (int_score=0).
	"""
	def __init__(self, *args, haplotype={}, N=1, interaction_func=max, **kwargs):
		super().__init__(*args, N=N, interaction_func=interaction_func, **kwargs)
		self.haplotype = {key: float(value) for key, value in haplotype.items()}

	def calc(self, sample_data):
		"""From Sharp2019: For haplotypes with an interaction the beta is taken from Table S3, without an interaction it is scored independently for each haplotype of the pair (Table S1)."""
		prs_score = super(Interaction, self).calc(sample_data=sample_data)
		int_score = self.calc_interaction(sample_data=sample_data)
		hap_score = 0 if int_score else self.calc_haplotype(sample_data=sample_data)
		return prs_score + int_score + hap_score

	def calc_haplotype(self, sample_data):
		"""Calculate the haplotype component from Sharp2019 TableS1. Should only be done when *no* interaction is present. See Figure S2."""
		hap_score = 0.0
		for allele, genotypes in sample_data.items():
			# For simplicity, assume genotype is a tuple of alleles (e.g., (0, 1) or (1, 1))
			if allele in self.haplotype:
				hap_score += self.haplotype.get(allele, 0) * sum(genotypes)
				logger.debug(f"calc_haplotype: {allele} found. Current Sum = {hap_score}")
		logger.info(f"calc_haplotype: Haplotype Sum = {hap_score}")
		return hap_score / self.N

	@classmethod
	def register_haplotype(cls, variant, haplotype, build=None):
		"""
		Little helper to register the weights in the haplotype tree
		
		variant (class): A variant from PGS
		haplotype (dict): The dict we're building
		"""
		from pksnp import Allele

		# TODO: This is a hack, and should be removed if we ever get the 'is_interaction' 'is_haplotype' stuff going.
		import re
		rsid = re.search("rs[0-9]+", variant.rsID)[0]

		haplotype[Allele(chromosome=variant.chr_name, position=variant.chr_position, rsid=rsid, allele=variant.effect_allele, build=build)] = float(variant.effect_weight)
