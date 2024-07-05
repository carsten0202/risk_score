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
	
	haplotype: Haplotype data from Table S1 in Sharp2019 et al. To be used *only* when no interaction is present (int_score=0).
	"""
	
	# Dictionary to translate Sharps haplotypes (Interaction)
	complex_alleles = {
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*03:0X;HLA-DQB1*03:02": "C:G",
		"HLA-DQA1*03:0X;HLA-DQB1*03:01_x_HLA-DQA1*01:0X;HLA-DQB1*05:01": "T:A",
		"HLA-DQA1*03:0X;HLA-DQB1*03:02_x_HLA-DQA1*03:0X;HLA-DQB1*03:02": "G/G:G/G", 
		"HLA-DQA1*03:0X;HLA-DQB1*03:01_x_HLA-DQA1*03:0X;HLA-DQB1*03:02": "T:G",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*03:0X;HLA-DQB1*03:01": "C:T",
		"HLA-DQA1*03:0X;HLA-DQB1*03:01_x_HLA-DQA1*02:01;HLA-DQB1*02:02": "T:T",
		"HLA-DQA1*04:01;HLA-DQB1*04:02_x_HLA-DQA1*02:01;HLA-DQB1*02:02": "T:T",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*01:0X;HLA-DQB1*05:01": "C:A",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*02:01;HLA-DQB1*02:02": "C:T",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*05:01;HLA-DQB1*02:01": "C/C:C/C",
		"HLA-DQA1*01:02;HLA-DQB1*06:02_x_HLA-DQA1*01:02;HLA-DQB1*06:02": "T/T:T/T",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*05:05;HLA-DQB1*03:01": "C:C",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*01:03;HLA-DQB1*06:03": "C:T",
		"HLA-DQA1*03:0X;HLA-DQB1*03:02_x_HLA-DQA1*01:03;HLA-DQB1*06:03": "G:T",
		"HLA-DQA1*03:0X;HLA-DQB1*03:01_x_HLA-DQA1*03:0X;HLA-DQB1*03:01": "T/T:T/T",
		"HLA-DQA1*03:02;HLA-DQB1*03:03_x_HLA-DQA1*03:0X;HLA-DQB1*03:02": "A:G",
		"HLA-DQA1*03:0X;HLA-DQB1*03:01_x_HLA-DQA1*05:05;HLA-DQB1*03:01": "T:C",
		"HLA-DQA1*05:01;HLA-DQB1*02:01_x_HLA-DQA1*01:02;HLA-DQB1*06:02": "C:T",
	}

	# Dictionary to translate Sharps haplotypes (No interaction)
	haplotype_alleles = {
		"HLA-DQA1*03:0X;HLA-DQB1*03:02": "G",
		"HLA-DQA1*01:02;HLA-DQB1*06:02": "T",
		"HLA-DQA1*05:01;HLA-DQB1*02:01": "C",
		"HLA-DQA1*02:01;HLA-DQB1*02:02": "T",
		"HLA-DQA1*05:05;HLA-DQB1*03:01": "C",
		"HLA-DQA1*01:0X;HLA-DQB1*05:01": "A",
		"HLA-DQA1*03:0X;HLA-DQB1*03:01": "T",
		"HLA-DQA1*01:03;HLA-DQB1*06:03": "T",
		"HLA-DQA1*02:01;HLA-DQB1*03:03": "G",
		"HLA-DQA1*04:01;HLA-DQB1*04:02": "T",
		"HLA-DQA1*01:0X;HLA-DQB1*05:03": "G",
		"HLA-DQA1*03:02;HLA-DQB1*03:03": "A",
		"HLA-DQA1*01:02;HLA-DQB1*06:09": "A",
		"HLA-DQA1*01:03;HLA-DQB1*06:01": "A",
	}

	def __init__(self, *args, haplotype={}, N=1, interaction_func=max, **kwargs):
		super().__init__(*args, N=N, interaction_func=interaction_func, **kwargs)
		self.haplotype = {key: float(value) for key, value in haplotype.items()}

	@property
	def rsids(self):
		"""Return the ids of the weighted alleles"""
		rsdict = super().rsids
		rsdict.update({allele.rsid:None for allele in self.haplotype})
		return rsdict

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
		logger.info(f"calc_haplotype: Sample={sample_data.sample}, Haplotype Sum = {hap_score}")
		return hap_score / self.N

	@classmethod
	def register_haplotype(cls, variant, haplotype, build=None):
		"""
		Little helper to register the weights in the haplotype tree
		
		variant (class): A variant from PGS
		haplotype (dict): The dict we're building
		"""
		from pksnp import Allele
		effect_allele = cls.haplotype_alleles.get(str(variant.effect_allele), str(variant.effect_allele))
		haplotype[Allele(chromosome=variant.chr_name, position=variant.chr_position, rsid=variant.rsID, allele=effect_allele, build=build)] = float(variant.effect_weight)
