#!/home/fls530/anaconda3/bin/python

###########################################################
#
# ---%%%  PKSNPS.py: Module to hold the SNP class  %%%---
#

import math
import re
import sys
from collections import OrderedDict, namedtuple, UserDict

import pklib.pkcsv as csv






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
	def __init__(self, CHROM=None, POS=None, ID=None, data=None, *args, **kwargs):
		super().__init__(data, *args, **kwargs)
		if not any([CHROM, POS, ID is None]):
			if re.search("^[0-9xXyYmM]+:[0-9]+", ID):
				(CHROM,POS,*_) = ID.split(":")
		self["CHROM"] = str(CHROM) if CHROM else None
		self["POS"]   = int(POS) if POS else None
		assert any([ID is not None, self["CHROM"] and self["POS"]]), "A Locus must have an ID or a Chromosome Position (Chrom + pos)."
		self["ID"]    = str(ID) if ID is not None else self.posid()
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

	def __setattr__(self, attr, value):
		super().__setattr__(attr, value)

	def posid(self):
		if all([self.CHROM, self.POS]):
			return self.CHROM + ":" + str(self.POS)
		return None



##################################################
#
# --%%  RUN: 'SNP' Class Definition  %%--

# Hmm... Not sure if 'genotype' is the right fit here... Maybe extend SNP into a pVCF class with genotype data?
class SNP(Locus):
	"""A simple class holding values for one biallelic SNP"""
	def __init__(self, *args, REF=None, ALT=None, INFO=None, genotype=None, **kwargs):
		super().__init__(*args, **kwargs)
		self["REF"]   = str(REF) if REF is not None else None
		self["ALT"]   = str(ALT) if ALT is not None else None
		self["INFO"]  = INFO if isinstance(INFO,dict) else dict() # Store INFO/TAG values similar to a VCF file
		self["genotype"] = genotype # This guy should be a list of namedtuple('genotype', 'subjectid bases dosage phased')

	def __eq__(self, other):
		# Check rsID == rsid?
		# What about build?
		same = [super().__eq__(other)]
		if isinstance(other, SNP):
			same.extend([other.REF in [self.REF, self.ALT]] if other.REF is not None else [])
			same.extend([other.ALT in [self.REF, self.ALT]] if other.ALT is not None else [])
		return all(same)

	def drop_genotype(self):
		self["genotype"] = None
		return self

	def genotype(self):
		"""Returns a dict of SubjectID:GenoType based on dosages from the self.genotype slot."""
		out = dict()
		for gt in self.get("genotype"):
			p=Dosage(ref=self.REF, alt=self.ALT, dosage=gt.dosage).p(gt.bases)
			out[gt.subjectid] = GenoType(ID=self.ID, CHROM=self.CHROM, POS=self.POS, genotype=gt.bases, phased=gt.phased, p=p)
		return out

	def getGenoType(self, dosage=None):
		"""Calculate one GenoType based on provided dosage."""
		if dosage is not None:
			p = Dosage(ref=self.REF, alt=self.ALT, dosage=dosage).p()
			return GenoType(ID=self.ID, CHROM=self.CHROM, POS=self.POS, genotype = self.REF * abs(self.mladd(dosage) - 2) + self.ALT * self.mladd(dosage), p=p)
		else:
			sys.exit("SNP Error: You must specify some kind of dosage-like score to convert to genotype.")

	@staticmethod
	def mladd(dosage):
		return round(float(dosage))




##################################################
#
# --%%  RUN: 'Allele' Class Definition  %%--

# So, here's the justification for the Allele class... we could make ie hashable!!!
class Allele(Locus):
	"""An allele is basically a locus with a base assigned. But it can also hold other things, like allelic risk scores"""
	def __init__(self, *args, allele="", dosage=1, **kwargs):
		super().__init__(*args, data={'allele':str(allele)}, **kwargs)
		self["dosage"] = float(dosage) # The allelic dosage; ie probability that the allele is present. Default: 1
		if self["CHROM"] and self['POS']: # Because posid is required if a position is given. Otherwise we are not hashable.
			self["ID"] = self.posid()
		assert self["allele"], "Allele must be a valid string of length > 0. For Deletions: Custom is to use the base just prior to the deletion."
		assert self["ID"], "Allele must have a valid ID of length > 0 that is hashable."

	def __contains__(self, other):
		sys.exit("__contains__ Not Implemented yet ;-)")
		return False

	def __eq__(self, other):
		if isinstance(other, Allele):
			return self.ID == other.ID and self.allele == other.allele
		return NotImplemented

	def __hash__(self):
		return hash((self.ID, self.allele))

	def __repr__(self):
		return self.ID

	def posid(self):
		out = super().posid()
		if out is not None:
			return out + ":" + self.allele
		return None




##################################################
#
# --%%  RUN: 'GenoType' Class Definition  %%--

class GenoType(Locus):
	"""Similar to an Allele, but no allele slot and the new slot 'genotype' is a tuple."""
	def __init__(self, *args, genotype=None, phased=False, p=None, **kwargs):
		try: geno = tuple(str(g) for g in genotype if str(g) not in [':','/','|','_'])
		except TypeError:
			sys.exit("GenoType Error: Provided genotype not iterable.")
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

	def getAlleles(self, useID=False):
		# We should probably do a real calc on dosage here...
		dosage = 2 if self.homozygous else 1
		if useID or None in (self.CHROM, self.POS):
			return [Allele(ID=self.ID, allele=allele, dosage=dosage) for allele in self.genotype]
		return [Allele(CHROM=self.CHROM, POS=self.POS, allele=allele, dosage=dosage) for allele in self.genotype]

	@property
	def homozygous(self):
		return all([gt == self.genotype[0] for gt in self.genotype])



#################################################
#
# --%%  Constructor Functions  %%--

def ReadGeno(genofobj, info):
	"""Input: File object to a geno file from SNPextractor.py; Return: A dict (subjectid:genotype_list)"""
	geno = OrderedDict()
	for subject in csv.DictReader(genofobj):
		subjectid = subject.pop("")
		geno[subjectid] = OrderedDict()
		for snpid,value in subject.items():
			dosage = sum([float(x) for x in re.split("[/|]+", value)])
			geno[subjectid][snpid] = info[snpid].getGenoType(dosage)
	return geno

# Should add this guy to SNP class since it returns a SNP object
def ReadInfo(infofobj):
	"""Input: File object to an info file from SNPextractor.py"""
	snps = OrderedDict()
	infoiter = csv.DictReader(infofobj)
	for info in infoiter:
		alt = re.sub("[\[\]]+","",info.get("ALT", ""))
		snps[info.get("ID")] = SNP(ID=info.get("ID"), CHROM=info.get("CHROM"), POS=info.get("POS"), REF=info.get("REF"), ALT=alt)
	return snps

# Should add this guy to SNP class since it returns a SNP object
def ReadVCF(vcfiter, drop_genotypes=True):
	"""Read SNP information from a VCF file using PyVCF. Because this is slow, sample/genotype information is not read by default, but can be enabled"""
	import vcf
	snps = OrderedDict()
	for record in vcfiter:
		if record.ALT == "P": # Pretty dirty hack to capture data which are not a locus; like haplotypes from snp2hla.
			snp = SNP(ID=record.ID, REF=record.REF, ALT=record.ALT, INFO=record.INFO, FORMAT=record.FORMAT)
		else:
			snp = SNP(ID=record.ID, CHROM=record.CHROM, POS=record.POS, REF=record.REF, ALT=record.ALT, INFO=record.INFO, FORMAT=record.FORMAT)
		if drop_genotypes is False:
# This should be part of a class method; not freeform like this... set_genotype, add_genotype... or someshit...
#    Maybe even call some class which extends SNP and not SNP itself...
			genotype = namedtuple('genotype', 'subjectid bases dosage phased')
			snp["genotype"] = []
			for gt in record.samples:
				dosage = gt["DS"] if "DS" in gt.data._fields else None
				bases = re.split("[/|]+", gt.gt_bases)
				snp["genotype"].append(genotype(subjectid=gt.sample, bases=bases, dosage=dosage, phased=gt.phased))
		snps[record.ID] = snp
	return snps



