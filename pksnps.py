#!/home/fls530/anaconda3/bin/python

###########################################################
#
# ---%%%  PKSNPS.py: Module to hold the SNP class  %%%---
#

import math
import pkcsv as csv
import re
import sys
from collections import OrderedDict, UserDict



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

#	def __setattr__(self, attr, value):
#		self[attr] = value



##################################################
#
# --%%  RUN: 'Allele' Class Definition  %%--

# So, here's the justification for the Allele class... we could make ie hashable!!!
class Allele(Locus):
	"""An allele is basically a locus with a base assigned. But it can also hold other things, like allelic risk scores"""
	def __init__(self, *args, allele="", **kwargs):
		super().__init__(*args, **kwargs)
		self["ALLELE"] = str(allele)
		self["ID"] += ":" + self["ALLELE"]
# We should probably make the above safer... and validate?

	def __eq__(self, other):
		if isinstance(other, Allele):
			return self.ID == other.ID
		elif isinstance(other, str):
			return self.ID == other
		return NotImplemented

	def __hash__(self):
		return hash(self.ID)



##################################################
#
# --%%  RUN: 'GenoType' Class Definition  %%--

class GenoType(Locus):
	"""Similar to an Allele, but the 'allele' slot is called 'genotype' and it's a list."""
	def __init__(self, *args, genotype=None, dosage=None, **kwargs):
		super().__init__(*args, **kwargs)
		try: iter(genotype)
		except TypeError:
			sys.exit("GenoType Error: Provided genotype not iterable.")
		else:
			self["GENOTYPE"] = [str(g) for g in genotype]
			self["ID"] += ":" + "".join(self["GENOTYPE"])

#	def __contains__(self, other):
#		return True

	def __eq__(self, other):
		if isinstance(other, GenoType):
			return self.ID == other.ID
		elif isinstance(other, str):
			return self.ID == other
		return NotImplemented

	def __hash__(self):
		return hash(self.ID)



##################################################
#
# --%%  RUN: 'SNP' Class Definition  %%--

class SNP(Locus):
	"""A simple class holding values for one biallelic SNP"""
	def __init__(self, *args, REF=None, ALT=None, info=None, **kwargs):
		super().__init__(*args, **kwargs)
		self["REF"]   = str(REF) if REF is not None else None
		self["ALT"]   = str(ALT) if ALT is not None else None
		self["INFO"]  = info if isinstance(info,dict) else dict() # Store INFO/TAG values similar to a VCF file

	def __eq__(self, other):
		same = [super().__eq__(other)]
		# Write here...
		# Check rsID == rsid
		if isinstance(other, SNP):
			same.extend([other.REF in [self.REF, self.ALT]] if other.REF is not None else [])
			same.extend([other.ALT in [self.REF, self.ALT]] if other.ALT is not None else [])
		# What about build?
		return all(same)

	def genotype(self, dosage=None):
		if dosage is not None:
			return GenoType(self.CHROM, self.POS, self.ID, genotype = self.REF * abs(self.mladd(dosage) -2) + self.ALT * self.mladd(dosage))
		else:
			sys.exit("SNP Error: You must specify some kind of dosage-like score to convert to genotype.")

	@staticmethod
	def mladd(dosage):
		return round(float(dosage))



#################################################
#
# --%%  Constructor Functions  %%--

def ReadInfo(infofile):
	snps = OrderedDict()
	with open(infofile) as f:
		infoiter = csv.DictReader(f)
		for info in infoiter:
			snps[info.get("ID")] = SNP(ID=info.get("ID"), CHROM=info.get("CHROM"), POS=info.get("POS"), REF=info.get("REF"), ALT=info.get("ALT"))
	return snps

def ReadRisk(riskfile):
	risks = OrderedDict()
	with open(riskfile) as f:
		riskiter = csv.DictReader(f)
		for risk in riskiter:
			chrom = risk.get("CHROM", risk.get("POSID","").split(":")[0])
			pos   = risk.get("POS", risk.get("POSID","").split(":")[1])
			rsid  = risk.get("RSID", risk.get("ID", risk.get("POSID", object())))
			beta  = risk.get("BETA", math.log(float(risk.get("ODDSRATIO", 1))))
			risks[rsid] = Allele(CHROM=chrom, POS=pos, ID=rsid, allele=risk.get("ALLELE"), BETA=beta)
	return risks



