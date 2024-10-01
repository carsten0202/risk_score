###########################################################
#
# ---%%%  sample_allele.py: From a VCF return the Sample Genotypes, one sample at a time  %%%---
#

from collections.abc import MutableMapping
import logging

from . import Allele, Call, Genotype

logger = logging.getLogger(__name__)

# Register the weighted alleles from Sharp2019 which relate to DQA1 & DQB1 haplotypes
SHARP_DQ_SNPs = [
	Allele(rsid='rs9275490',  allele='G'), Allele(rsid='rs17843689',  allele='T'),
	Allele(rsid='rs9273369',  allele='C'), Allele(rsid='rs17211699',  allele='T'),
	Allele(rsid='rs9469200',  allele='C'), Allele(rsid='rs10947332',  allele='A'),
	Allele(rsid='rs1281935',  allele='T'), Allele(rsid='rs62406889',  allele='T'),
	Allele(rsid='rs28746898', allele='G'), Allele(rsid='rs12527228',  allele='T'),
	Allele(rsid='rs1794265',  allele='G'), Allele(rsid='rs9405117',   allele='A'),
	Allele(rsid='rs16822632', allele='A'), Allele(rsid='rs117806464', allele='A'),
]

# There is no proper support for dosage scores in the script yet, but we could add it. Will likely need a dosage class:
# Will probably want to extend float on this one. Look here for some inspiration:
# https://stackoverflow.com/questions/25022079/extend-python-build-in-class-float

class HaplotypeError(Exception):
	pass

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

		alleles (dict): Allele:Call as list of genomic probabilities or dosage
		sample (str): Name of sample (optional)
		build (str): Genomic Build (optional)

		The alleles (dict) has an allele (rsid:base) as key and a list (of up to ploidy length). It can contain dosage
		scores or genotype probabilities. Same RSID can be given several times for different variants. Unobserved
		variants are typically not recorded, i.e. no '0' prob. entries.
		For genotype data it may look like this:
			{Allele(rs2476601, A): [1], Allele(rs2476601, G): [1], Allele(rs2111485, G): [1, 1]}
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

	@property
	def alleles(self):
		"""Returns the allele dict from _store."""
		return self._store

	def Sharp_DQ_markers(self, klitz:bool=False):
		"""Return SampleAlleles instance with only the snps used by Sharp2019 to infer HLA-DQA1-DQB1 haplotypes. Optional sorting according from Klitz2003."""
		# Ok, so we still fail with 'ignore' if the most frequent haplotype isn't weighted... See subject 58x1055
		sorted_alleles = {a:self[a] for a in SHARP_DQ_SNPs if self.get(a)}

		# Check number of occurrences for each allele are two-ish
		if sum(sorted_alleles.values()) > 2:
			if klitz:
				pairs = iter(sorted_alleles.items())
				(allele, call) = next(pairs)
				if call >= 2:
					sorted_alleles = {allele:call}
				else:
					(a2, c2) = next(pairs)
					sorted_alleles = {allele:call, a2:Call(0, GT=(0,1))}
					# Yes, by Klitz the second haplotype can only be heterozygous
				logger.debug(f" Klitz HLA indicators found = {sorted_alleles}")
			else:
				raise HaplotypeError(f"Spurious HLA Haplotype inferrence encountered for sample={self}")
		return SampleAlleles(alleles=sorted_alleles, build=self.build, sample=self.sample)
