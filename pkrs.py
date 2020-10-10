#!/home/fls530/anaconda3/envs/snakemake/bin/python

"""
Just like th 're' module; you pre-complie a GRS based on a dict of alt alleles and weights

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 

"""

import math
import pkcsv as csv
import pathlib

class RiskScore:
	"""A template object for holding the definition of a weight-based risk score"""
	def __init__(self, factors=None, alt_allele=dict()):
		self.allele = dict() # All of these should be dicts (OrderedDicts are acceptable); indexed by appropriate id (snpid?)
		self.beta   = dict()
		self.odds   = dict()
		if isinstance(factors, dict):
			self.N = len(factors)
			for fac in factors.values():
				self.append(fac)
		elif isinstance(factors, list):
			self.N = len(factors)
			for fac in factors:
				self.append(fac)
		else:
			raise SystemExit("Unable to process 'factor' argument. Please provide either a list-o-dicts or a dict-o-dicts")
		self.direct = {k:(v == alt_allele.get(k, v)) for k,v in self.allele.items()}

	def append(self, factor):
		facid = factor["SNPID"]
		self.allele[facid] = factor.get("ALLELE", None)
		self.odds[facid] = float(factor.get("ODDSRATIO", None))
		self.beta[facid] = float(factor.get("BETA", math.log(self.odds[facid])))


	def calc(self, snpdict):
		"""This function implements a simple risk score based on a weighted sum.
		   Subject should be a dict with str(id):float(allele)"""
		wsum = 0
		for snpid,allele in snpdict.items():
			wsum += self.beta.get(snpid,0) * abs(float(allele) - (0 if self.direct.get(snpid,True) else 2))
		return (wsum / self.N)


class oram2016(RiskScore):
	"""Oram et al 2016"""
	def weights():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt") as f:
			return list(csv.DictReader(f))

	def hlaweights():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights_hla.txt") as f:
			import collections
			dict_dicts = collections.defaultdict(lambda: collections.defaultdict(float))
			for line in csv.DictReader(f):
				odds = float(line.pop("ODDSRATIO",1))
				beta = float(line.pop("BETA", math.log(odds)))
				(id1,id2) = sorted(line.values())
				dict_dicts[id1][id2] = beta
			return dict_dicts

	def __init__(self, alt_allele=dict(), factors=weights(), hlafactors=hlaweights(), **kwargs):
		super().__init__(alt_allele=alt_allele, factors=factors, **kwargs)
		self.hla = hlafactors
		self.N = 2 * (self.N + 1)

	def calc_oramhla(self, snpdict):
		"""Calculate the HLA-part of Oram2016. Snpdict should have the format {str(id):float(allele)}"""
		wsum = 0
		for id1 in self.hla:
			for id2 in self.hla[id1]:
				snpid1,allele1 = id1.split("_")
				snpid2,allele2 = id2.split("_")
				if snpid1 in snpdict and int(allele1) == round(float(snpdict[snpid1])) and \
				   snpid2 in snpdict and int(allele2) == round(float(snpdict[snpid2])):
					wsum += self.hla[id1][id2]
		return (wsum / self.N)

	def calc(self, snpdict, **kwargs):
		wsum = super().calc(snpdict, **kwargs)
		return wsum + self.calc_oramhla(snpdict)


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


