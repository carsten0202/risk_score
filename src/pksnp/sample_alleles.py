###########################################################
#
# ---%%%  sample_allele.py: From a VCF return the Sample Genotypes, one sample at a time  %%%---
#

from collections.abc import MutableMapping
import logging

from . import Allele, Genotype

logger = logging.getLogger(__name__)


# There is no proper support for dosage scores in the script yet, but we could add it. Will likely need a dosage class:
# Will probably want to extend float on this one. Look here for some inspiration:
# https://stackoverflow.com/questions/25022079/extend-python-build-in-class-float

class SampleAlleles(MutableMapping):
	"""
	A dictionary that holds sample alleles from a VCF file

	Attributes:
	_store: A dictionary where keys are alleles and values are sample genotypes in tuples; e.g. (0,1) or for dosage e.g. (1.9).
	sample (str): Name of sample
	build (str): Genomic build
	"""
	def __init__(self, alleles=None, sample=None, build=None):
		"""
		Parses a VCF file and transposes the data to allow access by sample.

		alleles (dict): Allele:Genotype as list of genomic probabilities or dosage
		sample (str): Name of sample (optional)
		build (str): Genomic Build (optional)
		"""
		self.build = build
		self.sample = sample
		self._store = dict()
		self._store.update(alleles)

	def __contains__(self, other):
		"""Check if allele is registered in the sample and its genotype dosage is >0.5 (i.e. it is at least heterozygotic)"""
		if isinstance(other, Genotype):
			# OK, what to do here? We could split the Geno into Alleles and call this func again on those alleles (elegant...)
			# Won't work for hashing, but we are not using hashing on genotypes, only 'x in y'
			out = []
			cutoff = 0.5 + other.homozygote
			for allele in other.alleles:
				out.append(sum(self._store.get(allele, [0])) > cutoff) # Must be true for all alleles...
			return all(out)
		elif isinstance(other, Allele):
			return sum(self._store.get(other, [0])) > 0.5
		return NotImplemented

	def __delitem__(self, key):
		del self._store[key]

	def __getitem__(self, key):
		return self._store[key]

	def __repr__(self):
		return f"{self.sample}"

	def __setitem__(self, key, value):
		self._store[key] = value

	def __iter__(self):
		return iter(self._store)

	def __len__(self):
		return len(self._store)
