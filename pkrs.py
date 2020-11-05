#!/home/fls530/anaconda3/envs/snakemake/bin/python

"""
Just like the 're' module; you pre-complie a GRS based on a dict of alt alleles and weights

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 

"""

import collections
import math
import pathlib
import pkcsv as csv
import pksnps
import sys

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

# This guy should receive genotypes (Alleles?), or at least be able to...:
	def calc(self, gtdict):
		"""This function implements a simple risk score based on a weighted sum.
		   Subject should be a dict with str(id):float(genotype_score)"""
		wsum = 0
		for gtid,gt in gtdict.items():
			if isinstance(gt, pksnps.GenoType):
				for allele in gt.getAllele():
					wsum += self.beta.get(allele, 0) * allele.p
			else:
				sys.exit("ERROR: Looks like you called calc on something thats not a GenoType")
#				wsum += self.beta.get(gtid,0) * abs(float(gt) - (0 if self.direct.get(gtid,True) else 2))
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

#	def calcFromGeno(self, genofobj):
#		"""Input: File object to a geno file from SNPextractor.py"""
#		for subject in csv.DictReader(genofobj):
#			subjectid = subject.pop("")
#			geno = OrderedDict()
#			for snpid,value in subject.items():
#				dosage = sum([float(x) for x in re.split("[/|]+", value)])
#				geno[snpid] = GenoType(CHROM=info[snpid].CHROM, POS=info[snpid].POS, genotype=info[snpid].getGenoType(dosage), dosage=dosage)
#			print(subjectid + "\t" + self.calc(geno))



#################################################
#
# --%%  CLASS: MultiRiskScore  %%--

class MultiRiskScore(RiskScore):
	"""Calculate GRS based on MultiLocus Weights."""
	def __init__(self, snps, risks, multirisks, *args, **kwargs):
		super().__init__(snps=snps, risks=risks, *args, **kwargs)
		self.multi = self.ReadMultiRisk(multirisks)

	def calc(self, gtdict, **kwargs):
		"""INPUT: Gtdict => dict {str(id):GenoType}; RETURN: A risk score (float)"""
		wsum = super().calc(gtdict, **kwargs)
		for gt1 in self.multi:
			for gt2 in self.multi[gt1]:
				if gt1 in gtdict.values() and gt2 in gtdict.values():
					wsum += self.multi[gt1][gt2] / self.N
		return wsum

	@staticmethod
	def ReadMultiRisk(riskfobj):
		"""We should probably write this with arbitrary nesting... It's not completely finished and is therefore intentionally dirty"""
		dict_dicts = collections.defaultdict(dict)
		riskiter = csv.DictReader(riskfobj)
		for risk in riskiter:
			beta = float(risk.get("BETA", math.log(float(risk.get("ODDSRATIO")))))
			gt = []
			while True:
				i = len(gt) + 1
				chrom = risk.get("CHROM_" + str(i), risk.get("POSID_" + str(i),":").split(":")[0])
				pos   = risk.get("POS_" + str(i), risk.get("POSID_" + str(i),":").split(":")[1])
				try: gt.append(pksnps.GenoType(CHROM=chrom, POS=int(pos), genotype=risk.get("GENOTYPE_" + str(i), "").split(":")))
				except AttributeError as ae: 
					print("\n" + str(ae), file=sys.stderr)
					sys.exit("Read Error: Each line of '" + str(riskfobj.name) + "' must contain one weight and at least two recognizable positions and alleles.\n")
				if len(gt) >= 2:
					break
			dict_dicts[gt[0]][gt[1]] = beta
		return dict_dicts




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
#	sharp is in three parts (int / no-int) + others
		wsum = super().calc(gtdict, **kwargs)
		return wsum + self.calc_sharp(gtdict)

	def calc_sharp(self, gtdict):
		"""For haplotypes with an interaction the beta is taken from Table S3, without an interaction it is scored independently for each haplotype of the pair. This last part mustm be done here."""
		wsum = 0
		return wsum / self.N


