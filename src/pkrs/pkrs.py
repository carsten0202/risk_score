#!/home/fls530/anaconda3/bin/python

"""
Just like the 're' module; you pre-complie a GRS based on a dict of alt alleles and weights

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 

"""

import logging
import math
import re
import sys

import pksnp.pksnp as pksnp

from collections import OrderedDict, namedtuple, UserDict

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)



#################################################
#
# --%%  CLASS: RiskScore  %%--

class RiskScore:
	"""An algorithm template object for holding the definition of a weight-based risk score.
	   Input should be an iterable of SNPs and an iterable of Alleles (with BETA defined)"""
	def __init__(self, risks, N=None):
		self.beta   = dict()
		self.risks = self.ReadRisk(risks)
		self.N = len(self.risks) if N is None else float(N)
		assert self.N > 0, f"The denominator given with '-n' ('{self.N}') for the arithmetric mean must be >0."
		logger.debug(f"RiskScore: Setting N={self.N}")
		for risk in self.risks:
			self.beta[risk]   = float(risk.get("BETA",0))

	@property
	def alleles(self):
		"""Some kind of function to output the allelic information. i.e. the positions+bases, but not the weights."""
		return self._risks

	@property
	def risks(self):
		"""Alternate name for getting the risk alleles. Just for convenience."""
		return self._risks

	@property
	def risks_rsids(self):
		"""Get the rsids for the risks. Note that this guy is broken as the 'ID' slot may hold other types of ids.."""
		return [allele.ID for allele in self._risks]

	@risks.setter
	def risks(self, value):
		""""""
		self._risks = value

	def calc(self, alleles):
		"""This function implements a simple risk score based on a weighted sum.
			alleles: Subject dict with Allele:float(dosage)."""
		wsum = 0
		for allele, wgt in self.beta.items():
			wsum += wgt * alleles.get(allele, 0)
			logger.debug(f"Calc Aggregate: Allele = '{allele}'; weight = {wgt}, dosage = {alleles.get(allele,0)}; wsum = {wsum}")
		return wsum / self.N

	def calc_by_snp(self, snpiter):
		"""Calculate GRS incrementally by snp for all subjects. e.g. when reading from a vcf."""
		wsums = OrderedDict()
		for snp in snpiter:
			for subjectid, gt in snp.samples_iteritems():
				wsums[subjectid] = wsums.get(subjectid,0) + sum(self.beta.get(allele,0) * allele.dosage for allele in gt.alleles) / self.N
			logger.debug(f"{subjectid} [allele, dosage, risk] = {list([allele,allele.dosage,self.beta.get(allele,0)] for allele in gt.alleles)} wsum = {wsums[subjectid]}")
		return wsums

	@staticmethod
	def ReadRisk(riskiter):
		"""Input: An iterator giving something with a 'get' method; Return: List of Alleles"""
		risks = []
		logger.debug(f"ReadRisk: {type(riskiter)}")
		for risk in riskiter:
			logger.debug(f"ReadRisk: {risk}")
			rsid  = risk.get("RSID")
			chrom = risk.get("CHROM", re.split(":", risk.get("POSID",":"))[0])
			pos   = risk.get("POS", re.split(":|_", risk.get("POSID",":"))[1])
			beta  = risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1))))
			try: risks.append(pksnp.Allele(CHROM=str(chrom), POS=int(pos), allele=str(risk.get("ALLELE")), BETA=float(beta)))
			except AttributeError as ae: 
				print("\n" + str(ae), file=sys.stderr)
				sys.exit("Read Error: Each line of '" + str(riskiter.get("name","weights file")) + "' must contain at least weight value with a recognizable position and allele.\n")
			except ValueError as ve:
				print("\n" + str(ve), file=sys.stderr)
				sys.exit("Read Error: Looks like some weights data are not following the expected format.")
		return risks




#################################################
#
# --%%  CLASS: PGSCatalog  %%--

class PGSCatalog(RiskScore):
	"""Class to handle risk scores downloaded from the PGSCatalog."""
	def __init__(self, pgs, *args, risks=[], **kwargs):
		"""Like RiskScore; we just need to format the pgs risks correctly (pgs should be the return of pyclick.CSVFile() or
		   similar). Any predefined weights in 'risks' will be kept, mostly to provide a unified interface between the
		   RiskScore classes, but it does allow you to stack weights, if you want."""
		logger.debug(f"PGSCatalog: {type(pgs)}")
		for line in pgs:
			logger.debug(f"PGSCatalog: {line}")
			risk = {}
			risk['RSID']   = line.get('rsID')
			risk['CHROM']  = line.get('chr_name')
			risk['POS']    = line.get('chr_position')
			risk['ALLELE'] = line.get('effect_allele')
			risk['BETA']   = line.get('effect_weight')
			risks.append(risk)
		super().__init__(risks=risks, *args, **kwargs)

"""
PGS Format from: https://www.pgscatalog.org/downloads/#scoring_columns
rsID	dbSNP Accession ID (rsID)	Optional
chr_name	Location - Chromosome 	Required
chr_position	Location - Position within the Chromosome	Required
effect_allele	Effect Allele	Required
other_allele	Other allele(s)	Recommended
locus_name	Locus Name	Optional
is_haplotype
is_diplotype	FLAG: Haplotype or Diplotype	Optional
imputation_method	Imputation Method	Optional
variant_description	Variant Description	Optional
inclusion_criteria	Score Inclusion Criteria	Optional
effect_weight	Variant Weight	Required
is_interaction	FLAG: Interaction	Optional
is_dominant	FLAG: Dominant Inheritance Model	Optional
is_recessive	FLAG: Recessive Inheritance Model	Optional
dosage_0_weight	Effect weight with 0 copy of the effect allele	Optional
dosage_1_weight	Effect weight with 1 copy of the effect allele	Optional
dosage_2_weight
"""




#################################################
#
# --%%  CLASS: MultiRiskScore  %%--

class MultiRiskScore(RiskScore):
	"""Calculate GRS based on MultiLocus Weights."""
	def __init__(self, risks, multirisks, *args, **kwargs):
		super().__init__(risks=risks, *args, **kwargs)
		self.multi = self.ReadMultiRisk(multirisks)

	def calc(self, gtdict, **kwargs): # This guy still works on two levels; isn't nested like ReadMultiRisk now is.
		"""Calculate the Multilocus part of a GRS.
			gtdict => dict {str(snpid):GenoType};
			Return => A risk score (float)"""
		logger.debug(f"MultiRiskScore: Subject = {gtdict}")
		alleles = dict()
		for gt in gtdict.values():
			for allele in gt.alleles:
				alleles[allele] = allele.dosage
		wsum = super().calc(alleles, **kwargs)
		if mrs := self.nested_lookup(self.multi, gtdict.values()):
			logger.debug(f"MultiRiskScore: Found multirisk weight = {mrs}")
			wsum += mrs / self.N
		logger.debug(f"MultiRiskScore: Total score = {wsum}")
		return wsum

	@staticmethod
	def nested_lookup(nested_dict, subject): 
		"""Search the tree by looking for keys in the nested_dict
		nested_dict: 
		subject:     
		"""
		wsum = []
		if isinstance(nested_dict, dict):
			for haplo in nested_dict: # Pulling from nested ensures that the returned matching weight is the highest ranked (by fileorder); Also fast, only looping over existing keys.
				if haplo in subject:
					logger.debug(f"MultiRiskScore: Found allele '{haplo}' pointing to '{nested_dict[haplo]}'.")
#					return MultiRiskScore.nested_lookup(nested_dict[haplo], [gt for gt in subject if gt != haplo]) # Move down in nested structure. Exclude haplo from subject so it isn't counted again.
#					return MultiRiskScore.nested_lookup(nested_dict[haplo], subject) # Move down in nested structure.
					wsum.extend(sharp2019.nested_lookup(nested_dict[haplo], subject)) # Move down in nested structure.
		else:
			logger.debug(f"MultiRiskScore: Score Found = {nested_dict}")
			return [nested_dict] # Which should actually be the weight by now (a float)
		return wsum

	@staticmethod
	def ReadMultiRisk(risk_iter):
		"""INPUT: A file object to read from; RETURN: An arbitrarily nested dict with the required genotypes/haplotypes as keys and the weights as the bottommost values.
		   We should probably write this with arbitrary nesting... It's not completely finished and is therefore intentionally dirty"""
		def nested_read(risk, nested_dict=dict(), i=1):
			chrom = risk.get("CHROM_" + str(i), risk.get("POSID_" + str(i),":").split(":")[0])
			pos   = risk.get("POS_" + str(i), risk.get("POSID_" + str(i),":").split(":")[1])
			myid  = risk.get("ID_" + str(i))
			try:
				gtype = pksnp.GenoType(ID=myid, CHROM=chrom, POS=pos, genotype=risk.get("GENOTYPE_" + str(i), "").split(":"))
			except (AssertionError, AttributeError) as ae:
				beta = float(risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1)))))
				assert beta and isinstance(nested_dict, dict), "Each line of multilocus weights file must contain one weight and at least one recognizable allele or genotype."
				return beta
			nested_dict[gtype] = nested_read(risk, nested_dict.get(gtype, dict()), i+1)
			return nested_dict

		nested_dict = dict()
		for risk in risk_iter:
			logger.debug(f"ReadMultiRisk: Reading {risk}")
			if risk.get("GENOTYPE_1"):
				nested_dict = nested_read(risk, nested_dict)
		return nested_dict




#################################################
#
# --%%  CLASS: Oram2016  %%--

class oram2016(MultiRiskScore):
	"""Calculate GRS based on Oram et al 2016"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.N = 2 * (self.N + 1)




#################################################
#
# --%%  CLASS: Shapr2019  %%--

class sharp2019(MultiRiskScore):
	"""Calculate GRS based on Sharp et al 2019"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.N = 1

	def calc(self, gtdict, **kwargs):
		"""From Sharp2019: For haplotypes with an interaction the beta is taken from Table S3, without an interaction it is scored independently for each haplotype of the pair."""
		logger.debug(f"Sharp2019: Subject genotypes = {gtdict}")
		
		# Calculate the aggregate part
		alleles = dict()
		for gt in gtdict.values():
			for allele in gt.alleles:
				alleles[allele] = allele.dosage
		wsum = super(MultiRiskScore,MultiRiskScore).calc(self, alleles, **kwargs)

		# Only scan alleles which are present (dosage > 0.5)
		# Hack to fix some earlier problem with reading in data
		k = [k for k,v in alleles.items() if v > 0.5]
		subject = [allele for allele in k if allele.allele == 'P']
		logger.debug(f"Sharp2019: Scanning {len(subject)} types {subject}")

		# Calculate the multi-part
		if mrs := sum(self.nested_lookup(self.multi, subject)):
			logger.debug(f"Sharp2019: Found multirisk weights = {mrs}")
			wsum += mrs / self.N
		logger.debug(f"Sharp2019: Total score = {wsum}")
		return wsum

	@staticmethod
	def ReadMultiRisk(risk_iter):
		"""Returns: A nested dict with Genotypes and Alleles as keys and the matching weights as leaf values.
		risk_iter => An iterable containing tab-separated data on genotypes/alleles and their risk weights."""
		import itertools
		(risk_iter1, risk_iter2) = itertools.tee(risk_iter,2)
		nested_dict = super(sharp2019, sharp2019).ReadMultiRisk(risk_iter1)

		def nested_read(risk, nested_dict=dict(), i=1):
			chrom = risk.get("CHROM_" + str(i), risk.get("POSID_" + str(i),":").split(":")[0])
			pos   = risk.get("POS_" + str(i), risk.get("POSID_" + str(i),":").split(":")[1])
			myid  = risk.get("ID_" + str(i))
			try:
				gtype = pksnp.Allele(ID=myid, CHROM=chrom, POS=pos, allele=risk.get("ALLELE_" + str(i), ""))
			except (AssertionError, AttributeError) as ae:
				beta = float(risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1)))))
				assert beta and isinstance(nested_dict, dict), "Each line of multilocus weights file must contain one weight and at least one recognizable allele or genotype."
				return beta
			nested_dict[gtype] = nested_read(risk, nested_dict.get(gtype, dict()), i+1)
			return nested_dict

		for risk in risk_iter2:
			if risk.get("ALLELE_1"):
				nested_dict = nested_read(risk, nested_dict)
		return nested_dict



