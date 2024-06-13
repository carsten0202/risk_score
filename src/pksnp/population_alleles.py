###########################################################
#
# ---%%%  sample_allele.py: From a VCF return the Sample Genotypes, one sample at a time  %%%---
#

from collections import defaultdict
from collections.abc import MutableMapping
import logging
from pysam import VariantFile

from . import Allele, SampleAlleles

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
	def __init__(self, variantfile, build=None, filter_ids=None):
		"""
		Parses a VCF file and transposes the data to allow access by sample.

		Parameters:
		variantfile (VariantFile): VariantFile as returned by pysam.VariantFile(file_path)
		filter_ids (iterable): Only return entries where 'id in filter_ids'
		"""
		self._store = dict()
		self.build=build
		self.samples = variantfile.header.samples
		self.update({sample: defaultdict(list) for sample in self.samples})
		for record in variantfile:
			# rsid = translate_position_to_rsid(record.chrom, record.pos)
			# Auch!!! translate may get the wrong ID if there's overlapping features...
			if filter_ids and record.id not in filter_ids:
				next
			for sample in self:
				alleles = [record.ref] + list(record.alts)
				genotype = [alleles[nuc] for nuc in record.samples[sample]['GT']]
				assert len(genotype)==2, f"Sample {sample} is not biallelic!"
				for allele in genotype:
					self._store[sample][Allele(chromosome=record.chrom, position=record.pos, rsid=record.id, allele=allele)].append(1)

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
    dbsnp_file = '/projects/cbmr_shared/data/common_resources/dbsnp/b156/vcf/GCF_000001405.25.gz'

    dbsnp = VariantFile(dbsnp_file)
    chrom_names = dict(zip(range(1,25), get_contig_names(dbsnp_file)))

    # Fetch the record corresponding to the chromosome and position
    record = next(dbsnp.fetch(chrom_names[int(chromosome)], position - 1, position))
    
    if record:
        return record.id  # rsID
    else:
        return None
