#################################################
#
# --%%%  CLASS: Interaction  %%%--
#

import logging

logger = logging.getLogger(__name__)

from .riskscore import RiskScore

class Interaction(RiskScore):
	"""
	An algorithm template object for holding the definition of an interaction-based risk score.

	attributes:
	interaction: dict-o-dicts with weights as leafs
	interaction_func (func): function to handle when several weights are found in the tree
	"""
	def __init__(self, interaction=dict(), interaction_func=sum, *args, **kwargs):
		"""Forward some arguments to super-class
		interaction: A dict-o-dicts suitable for self.interaction
		interaction_func: Aggregate function for interaction_func
		"""
		super().__init__(*args, **kwargs)
		self.interaction = interaction
		self.interaction_func = interaction_func
		logger.debug(f" Weighted Interaction = {self.interaction}")

	def calc(self, sample_data):
		"""Overload calc to handle interactions"""
		prs_score = super().calc(sample_data=sample_data)
		int_score = self.calc_interaction(sample_data=sample_data)
		return prs_score + int_score

	def calc_interaction(self, sample_data):
		"""Calculate the interaction part of the score."""
		logger.debug(f"calc_interaction: Sample={sample_data.sample}, Scanning alleles={self.interaction}")
		try:
			int_score = self.interaction_func(self.traverse_interactions(sample_data))
			logger.info(f"calc_interaction: Sample={sample_data.sample}, Interaction Sum = {int_score}")
			return int_score / self.N
		except ValueError:
			logger.info(f"calc_interaction: Sample={sample_data.sample}, Interaction Sum = 0")
			return 0

	def traverse_interactions(self, sample_alleles):
		"""Function to traverse the tree and calculate PRS"""
		def recursive_traverse(node, score=[]):
			for key, child in node.items():
				if key == "weight":
					logger.debug(f"traverse: Found weight={child}")
					score.append(child)
				elif key in sample_alleles:
					logger.debug(f"traverse: Found {key}")
					score = recursive_traverse(child, score)
			return score

		logger.debug(f" Top Level={list(self.interaction.keys())}")
		return recursive_traverse(self.interaction)

	@classmethod
	def register_interaction(cls, variant, interaction):
		"""Little helper-function for FromPGS to register the weights in the interaction tree.
		
		variant (class): A variant from PGS
		interaction (dict): The dict we're building
		"""
		from pksnp import Allele, Genotype
		current_level = interaction
		branches = []

		rsids = cls.pat_rsid.findall(str(variant.rsID))
		genotypes =cls.pat_genotype.split(str(variant.effect_allele))
		for rsid,genotype in zip(rsids, genotypes):
			try:
				branches.append(Genotype(rsid=rsid, genotype=genotype))
			except IndexError:
				branches.append(Allele(rsid=rsid, allele=genotype))

		for branch in branches:
			if branch not in current_level:
				current_level[branch] = {}
			current_level = current_level[branch]
		current_level["weight"] = float(variant.effect_weight)

