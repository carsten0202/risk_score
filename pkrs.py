#!/home/fls530/anaconda3/envs/snakemake/bin/python

"""
Just like th 're' module; you pre-complie a GRS based on a dict of alt alleles and weights

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 

"""

import math
import pathlib
import pkcsv as csv
import pksnps

class RiskScore:
	"""An algorithm template object for holding the definition of a weight-based risk score.
	   Input should be an iterable of SNPs and an iterable of Alleles (with BETA defined)"""
	def __init__(self, snps, risks):
		self.beta   = dict()
		self.direct = dict()
		if isinstance(risks, dict):
			risks = risks.values()
		self.N = len(risks)
		if isinstance(snps, dict):
			snps = snps.values()
		self.snps   = dict(zip([s.ID for s in snps], snps))   # Dict of SNP instances
#		raise SystemExit("Unable to process 'factor' argument. Please provide either a list-o-dicts or a dict-o-dicts")
		for snp in snps:
			for risk in risks:
				if snp == risk:
					self.beta[snp.ID]   = risk.get("BETA",0)
					self.direct[snp.ID] = risk.get("ALLELE",object()) == snp.ALT

	def calc(self, gtdict):
		"""This function implements a simple risk score based on a weighted sum.
		   Subject should be a dict with str(id):float(genotype_score)"""
		wsum = 0
		for gtid,score in gtdict.items():
			wsum += self.beta.get(gtid,0) * abs(float(score) - (0 if self.direct.get(gtid,True) else 2))
		return wsum / self.N



class oram2016(RiskScore):
	"""Calculate GRS based on Oram et al 2016"""
	def __init__(self, snps, risks=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt", hlafactors=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights_hla.txt", **kwargs):
		super().__init__(snps=snps, risks=pksnps.ReadRisk(risks) if isinstance(risks, str) else risks, **kwargs)
		self.hla = self.hlaweights(hlafactors) # self.hla is a dict_dicts using 'GenoType's as keys.
		self.N = 2 * (self.N + 1)

	def calc_oramhla(self, gtdict):
		"""Calculate the HLA-part of Oram2016. Gtdict should have the format {str(id):float(allele)}"""
		wsum = 0
		gts  = [s.genotype(dosage=gtdict[i]) for i,s in self.snps.items()]
		gtsd = dict(zip([gt.ID for gt in gts],gts))
		for gt1 in self.hla:
			for gt2 in self.hla[gt1]:
				if gt1 in gtsd and gt2 in gtsd:
					wsum += self.hla[gt1][gt2]
		return (wsum / self.N)

	def calc(self, gtdict, **kwargs):
		wsum = super().calc(gtdict, **kwargs)
		return wsum + self.calc_oramhla(gtdict)

	def hlaweights(self, factors=None):
		import collections
		dict_dicts = collections.defaultdict(lambda: collections.defaultdict(float))
		if isinstance(factors, str):
			with open(factors) as f:
				for line in csv.DictReader(f):
					beta = float(line.get("BETA", math.log(float(line.get("ODDSRATIO")))))
					gt1 = pksnps.GenoType(CHROM=line.get("CHROM_1"), POS=line.get("POS_1"), hhh=line.get("ID_1"), genotype=line.get("GENOTYPE_1"))
					gt2 = pksnps.GenoType(CHROM=line.get("CHROM_2"), POS=line.get("POS_2"), hhh=line.get("ID_2"), genotype=line.get("GENOTYPE_2"))
					dict_dicts[gt1][gt2] = beta
		else:
			return NotImplemented
		return dict_dicts


class sharp2016(RiskScore):
	def weights():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.txt") as f:
			snptable = csv.DictReader(f)
			return list(csv.DictReader(f))

	def __init__(self, alt_allele=dict(), factors=weights(), **kwargs):
		super().__init__(alt_allele=alt_allele, factors=factors, **kwargs)
		self.N = 1

	def calc(self, subject):
#	sharp is in three parts (int / no-int) + others
		wsum = 0
		print(wsum)


