###########################################################
#
# ---%%%  sample_allele.py: From a VCF return the Sample Genotypes, one sample at a time  %%%---
#

from collections import defaultdict
from collections.abc import MutableMapping
import logging
from pysam import VariantFile

from . import Allele, Call, SampleAlleles

logger = logging.getLogger(__name__)


# There is no proper support for dosage scores in the script yet, but we could add it. Will likely need a dosage class:
# Will probably want to extend float on this one. Look here for some inspiration:
# https://stackoverflow.com/questions/25022079/extend-python-build-in-class-float

class PopulationAlleles(MutableMapping):
	"""
	A dictionary that holds sample alleles from a VCF file

	_store: A dictionary where keys are sample names and values are dictionaries of Allele objects and sample genotypes.
	samples: A list(str) with the names of samples
	"""
	def __init__(self, variantfile, build=None, filter_ids=None, format_field='GT'):
		"""
		Parses a VCF file and transposes the data to allow access by sample.

		Parameters:
		variantfile (VariantFile): VariantFile as returned by pysam.VariantFile(file_path)
		filter_ids (iterable): Only return entries where 'id in filter_ids'
		"""
		self._store = dict()
		self.variants = list()
		self.build = build
		self.samples = variantfile.header.samples
		self.update({sample: defaultdict(list) for sample in self.samples})

		# Read variants
		for record in variantfile:
			# if we add rsid translation, then we likely add it here...
			if filter_ids and record.id not in filter_ids:
				continue
			self.variants.append(record.id)
			alleles = [record.ref] + list(record.alts)
			for sample in self:
				for (i,a) in enumerate(alleles):
					self._store[sample][Allele(chromosome=record.chrom, position=record.pos, rsid=record.id, allele=a)]=Call(i, **{format_field:record.samples[sample][format_field]})

				# genotype = [alleles[nuc] for nuc in record.samples[sample]['GT']]
				# assert len(genotype)==2, f"ERROR: Sample={sample}, genotype={genotype} - is not biallelic!"
				# for allele in genotype:
				# 	self._store[sample][Allele(chromosome=record.chrom, position=record.pos, rsid=record.id, allele=allele)].append(1)

		# Validate variants
		logger.info(f" Read {len(self.variants)} variants from file '{variantfile.filename}")
		logger.debug(f" rsIDs found in input VCF = {self.variants}")
		overlap = [rsid in self.variants for rsid in filter_ids]
		if not all(overlap):
			no_missing = len(filter_ids)-sum(overlap)
			logger.warning(f" Of {len(filter_ids)} requested variants, {no_missing} were not found in input.")
			logger.info(f" Missing {no_missing} variants={[rsid for bool,rsid in zip(overlap, filter_ids) if bool is False]}")
			if no_missing > len(filter_ids) / 2:
				logger.error(f" More than half of the requested variants are missing. Likely something is very wrong!!!")
				logger.error(f" Check that your data files are the correct ones and that all contain rsIDs. rsIDs are necessary to correctly map variants.")
				import sys
				sys.exit("Exiting due to errors...")
		else:
			logger.info(f" All {len(filter_ids)} requested variants found in input data - Rejoice!")

	def __delitem__(self, key):
		del self._store[key]

	def __getitem__(self, key):
		return SampleAlleles(self._store[key], sample=key, build=self.build)

	def __setitem__(self, key, value):
		self._store[key] = value

	def __iter__(self):
		return iter(self._store)

	def __len__(self):
		return len(self._store)




# These two are kinda hack-ish... with hardlined paths, etc. Should have a cleaner and quicker approach through an indexed sql interface.
# They can also be wrong if(=WHEN!) there's overlapping features

def get_contig_names(vcf_file_path):
    """
    Reads and returns all contig names (chromosome names) from a VCF file.
    
    Parameters:
    vcf_file_path (str): The path to the VCF file.
    
    Returns:
    list: A list of contig names.
    """
    vcf = VariantFile(vcf_file_path)
    contig_names = list(vcf.header.contigs)
    return contig_names

def translate_position_to_rsid(chromosome, position, build='GRCh38'):
    # Open the dbSNP file, which should be indexed for fast access
	# TODO: Right now this is hardlined in (and GRCh37...), but it shouldn't be...
	# Auch!!! translate may get the wrong ID if there's overlapping features...
    dbsnp_file = '/projects/cbmr_shared/data/common_resources/dbsnp/b156/vcf/GCF_000001405.25.gz'

    dbsnp = VariantFile(dbsnp_file)
    chrom_names = dict(zip(range(1,25), get_contig_names(dbsnp_file)))

    # Fetch the record corresponding to the chromosome and position
    record = next(dbsnp.fetch(chrom_names[int(chromosome)], position - 1, position))
    
    if record:
        return record.id  # rsID
    else:
        return None
