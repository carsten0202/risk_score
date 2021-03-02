#!/home/fls530/miniconda3/envs/myscripts/bin/python

#  riskscore.py

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

import click
from collections import namedtuple
import logging
import pathlib
import re
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
 
import pkcsv as csv
import pkclick
import pksnps as snps
import pkrs as riskscore

Version = "1.0.1"

EPILOG = namedtuple('Options', ['fileformat','multiformat','legal'])(
fileformat = """

COLUMN-BASED DATAFILES:

Several options accept files containing data in tables/columns. Columns separators are auto-detected from the input and should work for tab-separated files, comma-seperated (csv) files, and space-separated files. 
The file must have one and only one headline as columns are identified by name. Any column name not recognized is ignored.

\b
Recognized names are:
ALLELE    - Identification of the weighted allele. Usually given as a nucleotide.
BETA      - The weight to be used. Used unmodified. Takes precedent over 'ODDSRATIO'.
CHROM     - Chromosome name or number.
ODDSRATIO - Ignored if column 'BETA' is given. Otherwise the natural logarithm of ODDSRATIO will be used as weights.
POS       - Chromosomal Position.
POSID     - An identifier of the type CHR:POS.

""",
legal = """
Written by Carsten Friis Rundsten <fls530@ku.dk>

LEGAL:

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
""",
multiformat = """
\b
Recognized names in multi-locus weights files:
ALLELE_#   - The weighted allele. Usually a nucleotide.
BETA       - The weight to be used. Takes precedent over 'ODDSRATIO'.
GENOTYPE_# - The weighted genotype. Usually given as nucleotide:nucleotide.
ID_#       - Locus identifier. Must match identifiers given in the subject data.
ODDSRATIO  - Ignored if column 'BETA' is given. Otherwise the natural logarithm of ODDSRATIO will be used as weights.

Replace '#' with a number of 1 or higher to indicate each locus group in a multi-locus weight.

""",
)

OPTION = namedtuple('Options', ['geno','info','log','vcf','weights','multiweights'])(
	geno = """Geno file of the type created by SNPextractor.""",
	info = """Info file of the type created by SNPextractor.""",
	log  = """Control logging. Valid levels: 'debug', 'info', 'warning', 'error', 'critical'.""",
	vcf  = """
Load VCF File using PyVCF VCFv4.0 and 4.1 parser for Python. Note that this option requires reading the entire VCF into
memory, which is not recommended for large VCF files. Use 'bcftools view --regions' or similar, to first reduce large
files in size.
""",
	weights = """Single locus risk weights file. See format description below on 'Column-Based Datafiles'.\n""",
	multiweights = """Multi-locus risk weights file. See format description below on 'Column-Based Datafiles'.\n"""
)

ScriptPath = str(pathlib.Path(__file__).resolve().parent.absolute())

# Notes and TODOs:

# TODO: Thees a couple of functions which could be cleaned at the bottom
# TODO: When reading from VCF we should include an option for selecting the field you want; eg GT or DS, etc.

# --%%  END: Perform Basic Setup  %%--
#
##################################################



##################################################
#
# --%%  RUN: Define Commands  %%--

class StdCommand(click.Command):
	def __init__(self, *args, epilog=None, **kwargs):
		super().__init__(*args, **kwargs)
		self.params.insert(0, click.Option(('-i','--info',), type=click.File(), help=OPTION.info))
		self.params.insert(0, click.Option(('-g','--geno',), type=click.File(), help=OPTION.geno))
		self.params.insert(0, click.Option(('--vcf',), type=pkclick.gzFile(mode='rb'), help=OPTION.vcf))
		self.epilog = EPILOG.fileformat + EPILOG.multiformat + EPILOG.legal

@click.group()
@click.version_option(version=Version)
@click.option('--log', default="warning", help=OPTION.log, show_default=True)
def main(log):
	"""Calculate the Genetic Risk Score (GRS) for a list of subjects based on predefined risk weights."""
	try:
		log_num = getattr(logging, log.upper())
	except AttributeError:
		raise ValueError(f"Invalid log level: '{log}'")
	logging.basicConfig(level=log_num)





@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-w','--weights', type=click.File(), default=None, help=OPTION.weights)
def aggregate(geno, info, vcf, weights):
	"""Calculate Aggregated Risk Score from user-provided weights."""
	(sbjgeno, snpinfo, weights) = process_args(geno, info, vcf, weights)
	grs = riskscore.RiskScore(snpinfo, risks=weights)
	for sbjid,gt in sbjgeno.items():
		print(sbjid + "\t" + str(grs.calc(gt)))

@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-m','--multilocus', type=click.File(), default=f"{ScriptPath}/oram2016.weights.multilocus.txt", show_default=True, help=OPTION.multiweights)
@click.option('-w','--weights', type=click.File(), default=f"{ScriptPath}/oram2016.weights.txt", show_default=True, help=OPTION.weights)
def oram2016(geno, info, vcf, weights, multilocus):
	"""Calculate Gene Risk Score based on Oram et al 2016.

REFERENCE:

\b
A type 1 diabetes genetic risk score can aid discrimination between type 1 and
type 2 diabetes in young adults.
RA Oram, K Patel, A Hill, B Shields, TJ McDonald, A Jones, AT Hattersley,
MN Weedon.
Diabetes care 39 (3), 337-344.
https://doi.org/10.2337/dc15-1111
"""
	(sbjgeno, snpinfo, weights) = process_args(geno, info, vcf, weights, multilocus)
	grs = riskscore.oram2016(snpinfo, risks=weights, multirisks=multilocus)
	for sbjid,gt in sbjgeno.items():
		print("{subject}\t{grs}".format(subject=sbjid, grs=round(grs.calc(gt),4)))

@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-m', '--multilocus', type=click.File(), default=f"{ScriptPath}/sharp2019.weights.multilocus.txt", help=OPTION.multiweights, show_default=True)
@click.option('-w', "--weights", type=click.File(), default=f"{ScriptPath}/sharp2019.weights.txt", help=OPTION.weights, show_default=True)
def sharp2019(geno, info, vcf, weights, multilocus):
	"""Calculate Gene Risk Score based on Sharp et al 2019.

REFERENCE:

\b
Development and standardization of an improved type 1 diabetes genetic risk
score for use in newborn screening and incident diagnosis.
SA Sharp, SS Rich, AR Wood, SE Jones, RN Beaumont, JW Harrison, DA Schneider,
JM Locke, JT, MN Weedon, WA Hagopian, RA Oram.
Diabetes Care 2019 Feb; 42(2): 200-207.
https://doi.org/10.2337/dc18-1785
"""
	(sbjgeno, snpinfo, weights) = process_args(geno, info, vcf, weights, multilocus)
	grs = riskscore.sharp2019(snpinfo, risks=weights, multirisks=multilocus)
	for sbjid,gt in sbjgeno.items():
		print("{subject}\t{grs}".format(subject=sbjid, grs=round(grs.calc(gt),4)))

@main.command(no_args_is_help=True, hidden=True)
@click.option('-m','--multilocus',  type=click.File(), default=f"{ScriptPath}/oram2016.weights.multilocus.txt", help=OPTION.multiweights, show_default=True)
@click.option('-w','--weights', type=click.File(), default=f"{ScriptPath}/oram2016.weights.txt", help=OPTION.weights, show_default=True)
def test(vcf, weights, multilocus):
	"""FOR TESTING PURPOSES ONLY; DO NOT USE!"""
	import vcf as vcf_
	vcfdata = snps.ReadVCF(vcf_.Reader(vcf, compressed=False), drop_genotypes=False)
	sbjgeno = transpose_geno([snp.genotype() for snp in vcfdata.values()])
	snpinfo = [snp.drop_genotype() for snp in vcfdata.values()]
	grs = riskscore.oram2016(snpinfo, risks=weights, multirisks=multilocus)
	for subjectid, gt in sbjgeno.items():
		print(subjectid + "\t" + str(grs.calc(gt)))

# --%%  END: Define Commands  %%--
#
##################################################




##################################################
#
# --%%  RUN: Subroutines  %%--

def process_args(geno=None, info=None, vcf=None, weights=None, *args):
	if vcf:
		import vcf as vcf_
		vcfdata = snps.ReadVCF(vcf_.Reader(vcf, compressed=False), drop_genotypes=False)
		sbjgeno = transpose_geno([snp.genotype() for snp in vcfdata.values()])
		snpinfo = [snp.drop_genotype() for snp in vcfdata.values()]
	elif (geno and info):
		snpinfo = snps.ReadInfo(info)
		sbjgeno = snps.ReadGeno(geno, snpinfo)
	else:
		sys.exit("ERROR: No input files specified, add either '--vcf' or both '--geno' and '--info' options. Use '--help' for help.")
	logging.info(f"Args: Reading weights from {weights.name}")
	weights = csv.DictReader(weights)
	logging.debug(f"{type(weights)}")
	return (sbjgeno, snpinfo, weights)

# Another quick'n'dirty function which should probably be done cleaner
def transpose_geno(genoiter):
	import collections
	dict_dicts = collections.defaultdict(dict)
	for geno_dict in genoiter:
		for sbjid, gt in geno_dict.items():
			dict_dicts[sbjid][gt.ID] = gt
	return dict_dicts

# --%%  END: Subroutines  %%--
#
##################################################


if __name__ == '__main__':
	main()


