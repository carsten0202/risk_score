#!/home/fls530/anaconda3/bin/python

"""
Just like the 're' module; you pre-complie a GRS based on a dict of alt alleles and weights

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 

"""

import collections
import math
import pathlib
import sys

import pkcsv as csv
import pksnps

#################################################
#
# --%%  CLASS: RiskScore  %%--

class RiskScore:
	"""An algorithm template object for holding the definition of a weight-based risk score.
	   Input should be an iterable of SNPs and an iterable of Alleles (with BETA defined)"""
	def __init__(self, snps, risks):
		self.beta   = dict()
		self.direct = dict()
		if isinstance(snps, dict):
			snps = snps.values()
		self.snps   = dict(zip([s.ID for s in snps], snps))   # Dict of SNP instances
		risks = self.ReadRisk(risks)
		self.N = len(risks)
		for snp in snps:
			for risk in risks:
				if snp == risk:
					self.beta[risk]   = float(risk.get("BETA",0))
					self.direct[risk] = risk.get("allele",object()) == snp.ALT

	def calc(self, gtdict):
		"""This function implements a simple risk score based on a weighted sum.
		   gtdict:	Subject dict with str(var_id):float(dosage)"""
		wsum = 0
		for gt in gtdict.values():
			if isinstance(gt, pksnps.GenoType):
				for allele in gt.getAlleles():
					wsum += self.beta.get(allele, 0) * allele.p
			else:
				sys.exit("ERROR: Looks like you called calc on something that's not a GenoType")
		return wsum / self.N

	@staticmethod
	def ReadRisk(riskfobj):
		"""Input: A file_obj (Or similar iterator); Return: List of Alleles"""
		risks = []
		riskiter = csv.DictReader(riskfobj)
		for risk in riskiter:
			chrom = risk.get("CHROM", risk.get("POSID",":").split(":")[0])
			pos   = risk.get("POS", risk.get("POSID",":").split(":")[1])
			beta  = risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1))))
			try: risks.append(pksnps.Allele(CHROM=str(chrom), POS=int(pos), allele=str(risk.get("ALLELE")), BETA=float(beta)))
			except AttributeError as ae: 
				print("\n" + str(ae), file=sys.stderr)
				sys.exit("Read Error: Each line of '" + str(riskfobj.name) + "' must contain at least weight value with a recognizable position and allele.\n")
		return risks




#################################################
#
# --%%  CLASS: MultiRiskScore  %%--

class MultiRiskScore(RiskScore):
	"""Calculate GRS based on MultiLocus Weights."""
	def __init__(self, snps, risks, multirisks, *args, **kwargs):
		super().__init__(snps=snps, risks=risks, *args, **kwargs)
		self.multi = self.ReadMultiRisk(csv.DictReader(multirisks))

	def calc(self, gtdict, **kwargs): # This guy still works on two levels; isn't nested like ReadMultiRisk now is.
		"""Calculate the Multilocus part of a GRS.
		   gtdict => dict {str(id):GenoType}; RETURN: A risk score (float)"""
		wsum = super().calc(gtdict, **kwargs)
		wsum += self.nested_lookup(self.multi, gtdict.values()) / self.N
		return wsum

	@staticmethod
	def nested_lookup(nested_dict, subject): 
		if isinstance(nested_dict, dict):
			for haplo in nested_dict: # Pulling from nested ensures that the returned matching weight is the highest ranked (by fileorder); Also fast, only looping over existing keys.
				if haplo in subject:
					return MultiRiskScore.nested_lookup(nested_dict[haplo], [gt for gt in subject if gt != haplo]) # Move down in nested structure. Exclude haplo from subject so it isn't counted again.
		else:
			return nested_dict # Which should actually be the weight by now (a float)
		return 0


	@staticmethod
	def ReadMultiRisk(risk_iter):
		"""INPUT: A file object to read from; RETURN: An arbitrarily nested dict with the required genotypes/haplotypes as keys and the weights as the bottommost values.
		   We should probably write this with arbitrary nesting... It's not completely finished and is therefore intentionally dirty"""
		def nested_read(risk, nested_dict=dict(), i=1):
			chrom = risk.get("CHROM_" + str(i), risk.get("POSID_" + str(i),":").split(":")[0])
			pos   = risk.get("POS_" + str(i), risk.get("POSID_" + str(i),":").split(":")[1])
			myid  = risk.get("ID_" + str(i))
			try:
				gtype = pksnps.GenoType(ID=myid, CHROM=chrom, POS=pos, genotype=risk.get("GENOTYPE_" + str(i), "").split(":"))
			except (AssertionError, AttributeError) as ae:
				beta = float(risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1)))))
				assert beta and isinstance(nested_dict, dict), "Each line of multilocus weights file must contain one weight and at least one recognizable allele or genotype."
				return beta
			nested_dict[gtype] = nested_read(risk, nested_dict.get(gtype, dict()), i+1)
			return nested_dict

		nested_dict = dict()
		for risk in risk_iter:
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
		wsum = super(MultiRiskScore,MultiRiskScore).calc(self, gtdict, **kwargs)
		wsum += sum(self.nested_lookup(self.multi, gtdict.values())[:2]) / self.N
		return wsum

	@staticmethod
	def nested_lookup(nested_dict, subject):
		if wsum := super(sharp2019,sharp2019).nested_lookup(nested_dict, subject):
			return [wsum]
		wsum = []
		import itertools
		subject_alleles = list(itertools.chain(*[gt.getAlleles(True) for gt in subject]))
		if isinstance(nested_dict, dict):
			for haplo in nested_dict: # Pulling from nested ensures that the returned matching weight is the highest ranked (by fileorder); Also fast, only looping over existing keys.
				if haplo in subject_alleles:
					wsum.extend(sharp2019.nested_lookup(nested_dict[haplo], subject)) # Move down in nested structure.
		else:
			return [nested_dict] # Which should actually be the weight by now (a float)
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
				gtype = pksnps.Allele(ID=myid, CHROM=chrom, POS=pos, allele=risk.get("ALLELE_" + str(i), ""))
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



