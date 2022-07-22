#!/home/fls530/anaconda3/bin/python

###########################################################
#
# ---%%%  PKSNPS.py: Module to hold the SNP class  %%%---
#

import logging
import math
import re
import sys
from collections import OrderedDict, namedtuple, UserDict

import pysam

import pklib.pkcsv as csv

logger = logging.getLogger(__name__)



# There is no proper support for dosage scores in the script yet, but we could add it. Will likely need a dosage class:
# Will probably want to extend float on this one. Look here for some inspiration:
# https://stackoverflow.com/questions/25022079/extend-python-build-in-class-float

class Dosage(object):
	"""Class for holding dosage scores."""
	def __init__(self, ref, alt, dosage):
		self.alt = str(alt) # This is the counted base
		self.ref = str(ref) # This is the Zero base
		self._dosage = float(dosage) if dosage is not None else None

	@property
	def dosage(self):
		"""Get the dosage."""
		return self._dosage

	def p(self, bases=None):
		"""I need to be able to convert a dosage score to an array of p-values. But that's not trivial..."""
		bases = bases if bases is not None else [self.ref, self.alt]
		decimals = self.dosage % 1 if self.dosage is not None else 0 # decimals=0 gives p=[1,1] signifying a fixed genotype
		if decimals > 0.5: # Here, decimals must be given to an ALT base
			p = [1.0, decimals] if bases[0] == self.ref else [decimals, 1.0]
		else: # Here, 1 - decimals must be given to a REF base
			p = [1.0, 1.0 - decimals] if bases[0] == self.alt else [1.0 - decimals, 1.0]
		return p





##################################################
#
# --%%  : 'Locus' Class Definition  %%--

class Locus(UserDict):
	"""A base class to hold one genomic position"""
	def __init__(self, CHROM=None, POS=None, ID=None, data={}, *args, **kwargs):
		super().__init__(data, *args, **kwargs)
		if not any([CHROM, POS, ID is None]):
			if re.search("^[0-9xXyYmM]+:[0-9]+", ID):
				(CHROM,POS,*_) = ID.split(":")
		self["CHROM"] = str(CHROM) if CHROM else None
		self["POS"]   = int(POS) if POS else None
		assert any([ID is not None, self["CHROM"] and self["POS"]]), "A Locus must have an ID or a Chromosome Position (Chrom + pos)."
		self["ID"]    = str(ID) if ID is not None else self.posid
		# What about build?

	def __eq__(self, other):
		"""A match is True if Chr+Pos and ID matches One pair can be missing; either chr+pos or ID, but not both."""
		if isinstance(other, Locus):
			same = [self.ID == other.ID]
			same.extend([self.CHROM == other.CHROM and self.POS == other.POS] if None not in [self.CHROM, other.CHROM, self.POS, other.POS] else [])
			return any(same)
		return NotImplemented

	def __getattr__(self, attr):
		try:
			return self[attr]
		except KeyError:
			raise AttributeError

#	def __setattr__(self, attr, value):
#		if attr == 'data':
#		        super().__setattr__(attr, value)
#		else:
#		        self.__setitem__(attr, value)

	@property
	def posid(self):
		if all([self.CHROM, self.POS]):
			return self.CHROM + ":" + str(self.POS)
		return None



##################################################
#
# --%%  RUN: 'SNP' Class Definition  %%--

class SNP(Locus):
	"""A simple class holding values for one biallelic SNP including optional sample information (e.g. a line of VCF)."""
	def __init__(self, *args, REF=None, ALT=None, INFO=None, samples=None, **kwargs):
		super().__init__(*args, **kwargs)
		self["REF"]   = str(REF) if REF is not None else None
		self["ALT"]   = str(ALT) if ALT is not None else None
		self["INFO"]  = INFO if isinstance(INFO,dict) else dict() # Store INFO/TAG values similar to a VCF file
		self.samples = samples

	def __eq__(self, other):
		# Check rsID == rsid?
		# What about build?
		same = [super().__eq__(other)]
		if isinstance(other, SNP):
			same.extend([other.REF in [self.REF, self.ALT]] if other.REF is not None else [])
			same.extend([other.ALT in [self.REF, self.ALT]] if other.ALT is not None else [])
		return all(same)

	@property
	def samples(self):
		return self.__samplerecord

	@samples.setter
	def samples(self, value):
		"""This guy should end up as a list of namedtuple('genotype', 'subjectid bases dosage phased')."""
		genotype = namedtuple('genotype', 'subjectid bases dosage phased')
		self.__samplerecord = []
		if isinstance(value, pysam.libcbcf.VariantRecordSamples):
			for acc, sample in value.iteritems():
				dosage = sample["DS"] if "DS" in sample.keys() else None
				self.__samplerecord.append(genotype(subjectid=acc, bases=sample.alleles, dosage=dosage, phased=sample.phased))
		elif value is not None:
			sys.exit("Not Implemented Yet!")

	def samples_iteritems(self):
		"""Generator for key, value pairs where key is Subject ID and Value is a GenoType object."""
		for sample in self.samples:
			p = [2 - sample.dosage, sample.dosage]
			gt = GenoType(ID=self.ID, CHROM=self.CHROM, POS=self.POS, genotype=f"{self.REF}{self.ALT}", phased=sample.phased, p=p)
			yield (sample.subjectid, gt)

	def getGenoType(self, dosage=None):
		"""Calculate one GenoType based on provided dosage."""
		if dosage is not None:
			p = Dosage(ref=self.REF, alt=self.ALT, dosage=dosage).p()
			return GenoType(ID=self.ID, CHROM=self.CHROM, POS=self.POS, genotype = self.REF * abs(self.mladd(dosage) - 2) + self.ALT * self.mladd(dosage), p=p)
		else:
			sys.exit("SNP Error: You must specify some kind of dosage-like score to convert to genotype.")

	@staticmethod
	def mladd(dosage):
		"""The Maximum Likelihood on the dosage."""
		return round(float(dosage))



class SNPiterator(object):
	"""An iterator which iterates over a list of SNPs (such as a VCF file) returning one SNP object at a time."""
	def __init__(self, variter, converter=lambda x: x, filterfun=None):
		"""This is all about preparing the super-important __variter iterator attribute."""
		"""filterids: collection of ids which convert with set(); converter: function that converts elements from variter to a SNP."""
		super().__init__()
		if filterfun is None:
			self.__variter = self._mygenerator(variter, converter=converter)
		else:
			self.__variter = self._mygenerator(filter(filterfun, variter), converter=converter)

	def __iter__(self):
		"""Return the iterator"""
		return self.__variter

	def __next__(self):
		"""Return the next element from the iterator."""
		return next(self.__variter)

	@staticmethod
	def _mygenerator(variter, converter):
		"""This generator parses input from the original variter through converter and returns a SNP object."""
		for count, record in enumerate(variter):
			if not count % 5000 and count:
				logger.info(f"SNPiterator: Reading variant data. {count} variants read.")
			yield converter(record)

	@classmethod
	def FromPySAM(cls, variter, filterids=None, contig=[None], start=[None], end=None, **kwargs):
		"""Constructor from a PySAM VariantFile object. filterid: only return SNP if record.id is in this list of strings."""
		logger.info(f"FromPySAM: Reading from file {variter.filename}")
		from itertools import chain
		end = start if end is None else end
		converter = lambda x: SNP(ID=x.id, CHROM=x.chrom, POS=x.pos, REF=x.ref, ALT="".join(x.alts), INFO=x.info, FORMAT=x.format, samples=x.samples)
		filterfun = None if filterids is None else lambda x: x.id in filterids
		iters = list()
		for c, s, e in zip(contig, start, end):
			iters.append(variter.fetch(contig=c, start=s, end=e, **kwargs))
		return cls(variter=chain(*iters), filterfun=filterfun, converter=converter)




##################################################
#
# --%%  RUN: 'Allele' Class Definition  %%--

# So, here's the justification for the Allele class... we could make ie hashable!!!
class Allele(Locus):
	"""An allele is basically a locus with a base assigned. But it can also hold other things, like allelic risk scores"""
	def __init__(self, *args, allele="", dosage=1, force_id=False, **kwargs):
		super().__init__(*args, data={'allele':str(allele)}, **kwargs)
		self["dosage"] = float(dosage) # The allelic dosage; ie probability that the allele is present. Default: 1
		if self["ID"] in ['.','']: # Because "false" ids cannot make us hashable.
			self["ID"] = self.posid
		assert self["allele"], "Allele must be a valid string of length > 0. For Deletions: Custom is to use the base just prior to the deletion."
		assert self["ID"], "Allele must have a valid ID of length > 0 that is hashable together with the allele."
		self.__hashkey = hash((self.ID, self.allele)) if (force_id or self.posid is None) else hash((self.posid, self.allele))

	def __contains__(self, other):
		sys.exit("__contains__ Not Implemented yet ;-)")
		return False

	def __eq__(self, other):
		if isinstance(other, Allele):
			return self.ID == other.ID and self.allele == other.allele
		return NotImplemented

	def __hash__(self):
		return self.__hashkey

	def __repr__(self):
		return f"{self.ID}:{self.allele}"




##################################################
#
# --%%  RUN: 'GenoType' Class Definition  %%--

# Should the GnoType not be a list/tuple of Alleles?
class GenoType(Locus):
	"""Similar to an Allele, but no allele slot and the new slot 'genotype' is a tuple."""
	def __init__(self, *args, genotype=None, phased=False, p=None, **kwargs):
		"""genotype: An iterable that returns the bases, e.g. a string; phased: bool; p: The p-value of the call, defaults to '1'."""
		try: geno = tuple(str(g) for g in genotype if str(g) not in [':','/','|','_'])
		except TypeError:
			sys.exit("GenoType Error: Provided genotype is not iterable.")
		super().__init__(*args, **kwargs)
		self["genotype"] = geno
		self["phased"] = phased
		self["p"] = p if p is not None else [1] * len(self["genotype"])

	def __contains__(self, other):
		sys.exit("__contains__ Not Implemented yet ;-)")
		return False

	def __eq__(self, other):
		if isinstance(other, GenoType):
			if self.phased and other.phased:
				return self.genotype == other.genotype and super().__eq__(other)
			return sorted(self.genotype) == sorted(other.genotype) and super().__eq__(other)
		return NotImplemented

	def __hash__(self): # This one isn't completely correct...
		if self.phased:
			return hash((self.ID, self.genotype))
		else:
			return hash((self.ID, tuple(sorted(self.genotype))))

	def __repr__(self):
		phased = "|" if self.phased else "/"
		return f"Genotype({self.ID}:{phased.join(self.genotype)},p={self.p})"

	@property
	def alleles(self):
		"""Return the alleles of the genotype."""
		return (Allele(CHROM=self.CHROM, ID=self.ID, POS=self.POS, allele=self.genotype[i], dosage=self.p[i]) for i in range(0,len(self.genotype)))

	@property
	def homozygous(self):
		return all([gt == self.genotype[0] for gt in self.genotype])

	@property
	def p(self):
		"""The probability of the noted alleles; essentially the allelic dosages of each allele in the genotype."""
		return self["p"]

	@p.setter
	def p(self, value):
		"""Set the p array."""
		self.__setitem__('p', value)

	def alleles_dosage(self, dosage):
		"""Return the two alleles from a diploid genotype based on a dosage value."""
		p = [2-dosage,dosage]
		return (Allele(CHROM=self.CHROM, POS=self.POS, ID=None, allele=self.genotype[i], dosage=p[i]) for i in (0,1))




#################################################
#
# --%%  CLASS: Cohort  %%--

# This is a bad name. Maybe VariantCalls? Other ideas?
class Cohort(object):
	"""An Umbrella Class for holding a collection of SNPs from a cohort with or without Genotype data. Essentially all info from a Variant File."""
	def __init__(self, snpinfo=None, sbjgeno=None, *args, **kwargs):
		"""Init cohort object. snpinfo: [SNPs] and sbjgeno: {sbjid:{snpid:GenoType}} are snp-centric data."""
		super().__init__(*args, **kwargs)
		self._snpinfo   = snpinfo
		self._genotypes = sbjgeno

	@property
	def genotypes(self):
		"""Get the genotype information."""
		return self._genotypes

	@staticmethod
	def _transpose(genoiter):
		"""Build the genotype matrix (SNP-centric GenoTypes -> {sbjid:{snpid:GenoType}}) (often needed when reading from a VCF file)."""
		import collections
		dict_dicts = collections.defaultdict(dict)
		for geno_dict in genoiter:
			for sbjid, gt in geno_dict.items():
				dict_dicts[sbjid][gt.ID] = gt
		return dict_dicts

	@classmethod
	def FromPySAM(cls, variter, filterids=None):
		"""Constructor from a PySAM VariantFile object. filterid: only return SNP if record.id is in this list of strings."""
		snps = OrderedDict()
		snpinfo = []
		for count, record in enumerate(variter.fetch()):
			if filterids and record.id not in filterids:
				next
			if not count % 1000 and count:
				logger.info(f"FromPySAM: Reading variant data. {count} variants read.")
			if "P" in record.alts: # Pretty dirty hack to capture data which are not a locus; like haplotypes from snp2hla.
				snp = SNP(ID=record.id, REF=record.ref, ALT=record.alts, INFO=record.info, FORMAT=record.format)
			else:
				snp = SNP(ID=record.id, CHROM=record.chrom, POS=record.pos, REF=record.ref, ALT=record.alts, INFO=record.info, FORMAT=record.format)
			snpinfo.append(snp)
			snp.samples = record.samples
			snps[record.id] = snp
		sbjgeno = cls._transpose([dict(snp.samples_iteritems()) for snp in snps.values()])
		logger.info(f"FromPySAM: Finished reading. {count} variants read.")
		return cls(snpinfo, sbjgeno)




#################################################
#
# --%%  Constructor & Other Support Functions  %%--

def read_alleles(genofobj, infofobj):
	"""Open and prepare data from a geno/info file object pair using the classes in pksnp.
		Returns => Tuple with subjectid and the alleles, an OrderedDict with Allele:dosage."""
	proto_genotypes = dict()
	infoiter = csv.DictReader(infofobj)
	for info in infoiter:
		alt   = re.sub("[\[\]]+","",info.get("ALT", ""))
		snpid = info.get("ID")
		proto_genotypes[snpid] = GenoType(CHROM=info.get("CHROM"), POS=info.get("POS"), genotype=(info.get("REF"), alt))
	for subject in csv.DictReader(genofobj):
		alleles   = OrderedDict()
		subjectid = subject.pop("")
		for snpid, value in subject.items():
			dosage = sum([float(x) for x in re.split("[/|]+", value)])
			for allele in proto_genotypes[snpid].alleles_dosage(dosage):
				alleles[allele] = allele.dosage
		logger.debug(f"read_geno: Subject '{subjectid}' alleles: {alleles}")
		yield (subjectid, alleles)

def read_genotypes(genofobj, infofobj):
	"""Open and prepare data from a geno/info file object pair using the classes in pksnp.
		Returns => Tuple with subjectid and the genotypes, an OrderedDict with GenoType:dosage."""
	proto_genotypes = dict()
	infoiter = csv.DictReader(infofobj)
	for info in infoiter:
		alt   = re.sub("[\[\]]+","",info.get("ALT", ""))
		snpid = info.get("ID")
		proto_genotypes[snpid] = GenoType(CHROM=info.get("CHROM"), POS=info.get("POS"), genotype=(info.get("REF"), alt))
	for subject in csv.DictReader(genofobj):
		subjectgt = proto_genotypes
		subjectid = subject.pop("")
		for snpid, value in subject.items():
			dosage = sum([float(x) for x in re.split("[/|]+", value)])
			subjectgt[snpid].p = [2-dosage, dosage]
		logger.debug(f"read_genotypes: Subject '{subjectid}' genotypes: {subjectgt}")
		yield (subjectid, subjectgt)

def read_pysam(variter, *args, **kwargs):
	"""This guy is different, because we need to iterate over subjects and not snps here. For Oram & Sharp."""
	snpiter = SNPiterator.FromPySAM(variter, *args, **kwargs)
	sbjgeno = __transpose([dict(snp.samples_iteritems()) for snp in snpiter])
	for subjectid, subjectgt in sbjgeno.items():
		yield (subjectid, subjectgt)

def __transpose(genoiter):
	"""Transpose the genotype matrix (list of {sbjid:GenoType} -> {sbjid:{gtid:GenoType}}) (often needed when reading from a VCF file)."""
	import collections
	dict_dicts = collections.defaultdict(dict)
	for geno_dict in genoiter:
		for sbjid, gt in geno_dict.items():
			dict_dicts[sbjid][gt.ID] = gt
	return dict_dicts




# This guy needs cleaning...
def ReadPySAM(variter, drop_genotypes=True):
	"""Read SNP information from a BCF/VCF file using PySAM. Sample/genotype information is not read by default, but can be enabled"""
	snps = OrderedDict()
	for record in variter.fetch():
		if "P" in record.alts: # Pretty dirty hack to capture data which are not a locus; like haplotypes from snp2hla.
			snp = SNP(ID=record.id, REF=record.ref, ALT=record.alts, INFO=record.info, FORMAT=record.format)
		else:
			snp = SNP(ID=record.id, CHROM=record.chrom, POS=record.pos, REF=record.ref, ALT=record.alts, INFO=record.info, FORMAT=record.format)
		if drop_genotypes is False:
# This should be part of a class method; not freeform like this... set_genotype, add_genotype... or someshit...
#    Maybe even call some class which extends SNP and not SNP itself...
			genotype = namedtuple('genotype', 'subjectid bases dosage phased')
			snp["genotype"] = []
			for acc, sample in record.samples.iteritems():
				dosage = sample["DS"] if "DS" in sample.keys() else None
				bases = sample.alleles
				snp["genotype"].append(genotype(subjectid=acc, bases=bases, dosage=dosage, phased=sample.phased))
		snps[record.id] = snp
	return snps


