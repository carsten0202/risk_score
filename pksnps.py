#!/home/fls530/anaconda3/bin/python

###########################################################
#
# ---%%%  PKSNPS.py: Module to hold the SNP class  %%%---
#

import math
import pkcsv as csv
import re
import sys
from collections import OrderedDict, namedtuple, UserDict



##################################################
#
# --%%  : 'Locus' Class Definition  %%--

class Locus(UserDict):
	"""A base class to hold one genomic position"""
	def __init__(self, CHROM=None, POS=None, ID=None, *args, **kwargs):
		super().__init__(*args, **kwargs)
		if (CHROM is None or POS is None):
			if ID is None:
				sys.exit("Locus error. Both genomic position and id cannot be 'None'. Please specify at least one of them.")
			elif re.search("^[0-9xXyYmM]+:[0-9]+", ID):
				(CHROM,POS,*_) = ID.split(":")
		self["CHROM"] = str(CHROM) if CHROM is not None else None
		self["POS"]   = int(POS) if POS is not None else None
		self["ID"]    = str(ID) if ID is not None else self["CHROM"] + ":" + str(self["POS"])
		# What about build?
		if not any([self["ID"], self["CHROM"] and self["POS"]]):
			sys.exit("LOCUS Error: You must specify at least an ID or a Chromosome Position (Chrom + pos).")

	def __eq__(self, other):
		"""A match is True if either Chr+Pos and ID matchs. One pair can be missing; either chr+pos or ID, but not both."""
		if isinstance(other, Locus):
			same = [self.ID == other.ID]
			same.extend([self.CHROM == other.CHROM and self.POS == other.POS] if None not in [self.CHROM, other.CHROM, self.POS, other.POS] else [])
			return any(same)
		print("Unsupported type: " + str(type(other)))
		return NotImplemented

	def __getattr__(self, attr):
		try:
			return self[attr]
		except KeyError:
			raise AttributeError

	def __setattr__(self, attr, value):
		super().__setattr__(attr, value)



##################################################
#
# --%%  RUN: 'SNP' Class Definition  %%--

class SNP(Locus):
	"""A simple class holding values for one biallelic SNP"""
	def __init__(self, *args, REF=None, ALT=None, INFO=None, genotype=None, **kwargs):
		super().__init__(*args, **kwargs)
		self["REF"]   = str(REF) if REF is not None else None
		self["ALT"]   = str(ALT) if ALT is not None else None
		self["INFO"]  = INFO if isinstance(INFO,dict) else dict() # Store INFO/TAG values similar to a VCF file
		self["genotype"] = genotype

	def __eq__(self, other):
		same = [super().__eq__(other)]
		# Write here...
		# Check rsID == rsid?
		if isinstance(other, SNP):
			same.extend([other.REF in [self.REF, self.ALT]] if other.REF is not None else [])
			same.extend([other.ALT in [self.REF, self.ALT]] if other.ALT is not None else [])
		# What about build?
		return all(same)

	def drop_genotype(self):
		self["genotype"] = None
		return self

	def genotype(self):
		out = OrderedDict()
		for gt in self.get("genotype"):
			out[gt.subjectid] = GenoType(CHROM=self.CHROM, POS=self.POS, genotype=gt.bases, dosage=gt.dosage)
		return out

	def getGenoType(self, dosage=None):
		if dosage is not None:
			# Should probably do a more intelligent p-calc here...
			return GenoType(CHROM=self.CHROM, POS=self.POS, genotype = self.REF * abs(self.mladd(dosage) - 2) + self.ALT * self.mladd(dosage))
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
	def __init__(self, *args, allele="", p=1, **kwargs):
		super().__init__(*args, **kwargs)
		self["allele"] = str(allele)
		self["ID"] += ":" + self["allele"]
		self["p"] = float(p) # The probability that the allele is present. Default: 1
		if not self["allele"]:
			raise AttributeError("ALLELE Error: Allele must be a valid string of length > 0. Indel Hint: Custom is to use the base just prior to the deletion ;-).")

	def __contains__(self, other):
		sys.exit("__contains__ Not Implemented yet ;-)")
		return False

	def __eq__(self, other):
		if isinstance(other, Allele):
			return self.ID == other.ID
		elif isinstance(other, str):
			return self.ID == other
		return NotImplemented

	def __hash__(self):
		return hash(self.ID)

	def getAllele(self):
		return [self]



##################################################
#
# --%%  RUN: 'GenoType' Class Definition  %%--

class GenoType(Allele):
	"""Similar to an Allele, but the 'allele' slot is colon-separated, and the new slot 'genotype' is a list."""
	def __init__(self, *args, allele=None, genotype=None, dosage=None, **kwargs):
		try: geno = [str(g) for g in genotype]
		except TypeError:
			sys.exit("GenoType Error: Provided genotype not iterable.")
		else:
			super().__init__(*args, allele=":".join(geno), **kwargs)
			self["genotype"] = geno
			self["dosage"] = dosage # This may not be the best way; dosages are hard to intrepret here...

	def __contains__(self, other):
		sys.exit("__contains__ Not Implemented yet ;-)")
		return False

	def __eq__(self, other):
		if isinstance(other, GenoType):
			return all(allele in self.getAllele() for allele in other.getAllele()) and all(allele in other.getAllele() for allele in self.getAllele())
		return NotImplemented

	def __hash__(self):
		return super().__hash__()

	def getAllele(self):
		# We should probably do a real calc on p here...
		return [Allele(CHROM=self.CHROM,POS=self.POS,allele=allele,p=1) for allele in self.genotype]



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

def ReadInfo(infofobj):
	"""Input: File object to an info file from SNPextractor.py"""
	snps = OrderedDict()
	infoiter = csv.DictReader(infofobj)
	for info in infoiter:
		snps[info.get("ID")] = SNP(ID=info.get("ID"), CHROM=info.get("CHROM"), POS=info.get("POS"), REF=info.get("REF"), ALT=info.get("ALT"))
	return snps

def ReadVCF(vcfiter, drop_genotypes=True):
	"""Read SNP information from a VCF file using PyVCF. Because this is slow, sample/genotype information is not read by default, but can be enabled"""
	import vcf
	snps = OrderedDict()
	for record in vcfiter:
		snp = SNP(ID=record.ID, CHROM=record.CHROM, POS=record.POS, REF=record.REF, ALT=record.ALT, INFO=record.INFO, FORMAT=record.FORMAT)
		if drop_genotypes is False:
# This should be part of a class method; not freeform like this... set_genotype, add_genotype... or someshit...
			genotype = namedtuple('genotype', 'subjectid bases dosage phased')
			snp["genotype"] = []
			for gt in record.samples:
				dosage = gt["DS"] if "DS" in gt.data._fields else None
				bases = re.split("[/|]+", gt.gt_bases)
				snp["genotype"].append(genotype(subjectid=gt.sample, bases=bases, dosage=dosage, phased=gt.phased))
		snps[record.ID] = snp
	return snps



