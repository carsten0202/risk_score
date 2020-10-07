#!/home/fls530/anaconda3/envs/snakemake/bin/python

import math
import mycsv as csv
import pathlib

class RiskScore:
	"""A data object for holding the definition of a GRS score??"""
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


	def calc(self, subject):
		"""This function implements a simple risk score based on a weighted sum"""
		wsum = 0
		for snpid,dosage in subject.items():
			wsum += self.beta.get(snpid,0) * abs(float(dosage) - (0 if self.direct.get(snpid,True) else 2))
		return (wsum / self.N)


class oram2016(RiskScore):
	"""Oram et al 2016"""
	def weights():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt") as f:
			return list(csv.DictReader(f))

	def interactions():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.interactions.txt") as f:
			out = {}
			for line in csv.DictReader(f):
				beta = line["ODDSRATIO"]
				out["ID"] = line["ODDSRATIO"]
			return out

	def __init__(self, alt_allele=dict(), factors=weights(), interactions=interactions(), **kwargs):
		super().__init__(alt_allele=alt_allele, factors=factors, **kwargs)
		print(interactions)
		self.N = 2 * (self.N + 1)

	def calc_oramhla(arg):
		return 0

	def calc(self, subject):
		wsum = super().calc(subject)
		print(wsum)
		print(wsum + self.calc_oramhla())
		return wsum


class sharp2016(RiskScore):
	def weights():
		with open(str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.txt") as f:
			snptable = csv.DictReader(f)
			return list(csv.DictReader(f))

	def __init__(self, alt_allele=dict(), factors=weights(), **kwargs):
		super().__init__(alt_allele=alt_allele, factors=factors, **kwargs)
		self.N = 1

	def calc(self, subject):
		wsum = super().calc(subject)
		print(wsum)

#	sharp is in three parts (int / no-int) + others

#	for allele in alleles:
#		if subject has allele:
#			gsr += weight

# Just like with regular expressions; we pre-complie a GRS based on some weights and a function

# YES! Use regex as template:

# Pattern are the weights
# Compile function is the actual algorithm
# Results are obtained by applying (ie 'matching') 


