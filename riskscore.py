#!/home/fls530/anaconda3/bin/python

#
#  riskscore.py
#
#  Copyright 2019 Carsten Friis Rundsten <fls530@ku.dk>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#


# Notes and TODOs:

# TODO: Thees a couple of functions which could be cleaned at the bottom
# TODO: There's still some mess in how we match weights => Allele, GenoType, dosage 
# TODO: Here's a funny bug. If the info file is a csv, then the hla part of oram fails...


##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

Version = "0.6 (Development Version)"

import click
import pathlib
import pkcsv as csv
import pkclick
import pksnps as snps
import pkrs as riskscore
import re
import sys
 
# --%%  END: Perform Basic Setup  %%--
#
##################################################

EPILOG_FileFormat = """
\b
Column-Based Datafiles:

Several options accept files containing data in tables/columns. Columns separators are auto-detected from the input and should work for tab-separated files, comma-seperated (csv) files, and space separated files. 
The file must have one and only one headline as columns are identified by name. Any column name not recognized is ignored.

\b
Recognized names are:
ALLELE    - Identification of the weighted allele. Usually given as a nucleotide. 
BETA      - The weight to be used. Used unmodified. Takes precedent over 'ODDSRATIO'.
CHROM     - Chromosome name or number.
ODDSRATIO - Ignored if column 'BETA' is given. Otherwise the natural logarithm of ODDSRATIO will be used as weights.
POS       - Chromosomal Position.
POSID     - An identifier of the type CHR:POS.

"""

OPTION_geno = """
Geno file of the type created by SNPextractor
"""

OPTION_info = """
Info file of the type created by SNPextractor
"""

OPTION_vcf = """
Load VCF File using PyVCF VCFv4.0 and 4.1 parser for Python. Note that this option requires reading the entire VCF into memory, which is not
recommended for large VCF files. Use 'bcftools view --regions' or similar, to first reduce large files in size.
"""

OPTION_weights = """
Single locus risk weights file. For accepted file formats, see description on 'Column-Based Datafiles'
"""

##################################################
#
# --%%  RUN: Commands  %%--

class StdCommand(click.Command):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.params.insert(0, click.Option(('-i','--info',), type=click.File(), help=OPTION_info))
		self.params.insert(0, click.Option(('-g','--geno',), type=click.File(), help=OPTION_geno))
		self.epilog = EPILOG_FileFormat

@click.group()
@click.version_option(version=Version)
def main():
	"""THIS SCRIPT IS STILL EXPERIMENTAL; USE WITH CAUTION"""
	pass

@main.command(cls=StdCommand)
@click.option('-w','--weights', type=click.File(), default=None,
              help="Single locus risk weights file.")
def aggregate(geno, info, weights):
	"""Calculate Aggregate Gene Risk Score from user-provided weights."""
	check_args(geno, info, weights)
	snpinfo = snps.ReadInfo(info)
	sbjgeno = snps.ReadGeno(geno, snpinfo)
	grs = riskscore.RiskScore(snpinfo, risks=weights)
	for sbjid,gt in sbjgeno.items():
		print(sbjid + "\t" + str(grs.calc(gt)))

@main.command(cls=StdCommand)
@click.option('-m','--multilocus',  type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.multilocus.txt", 
              help="Multilocus risk weights.")
@click.option('-w','--weights', type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt", help=OPTION_weights)
def oram2016(geno, info, weights, multilocus):
	"""Calculate Gene Risk Score based on Oram et al 2016."""
	check_args(geno, info, weights, multilocus)
	snpinfo = snps.ReadInfo(info)
	sbjgeno = snps.ReadGeno(geno, snpinfo)
	grs = riskscore.oram2016(snpinfo, risks=weights, multirisks=multilocus)
	for sbjid,gt in sbjgeno.items():
		print(sbjid + "\t" + str(grs.calc(gt)))

@main.command(cls=StdCommand)
@click.option('-m','--multilocus',  type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.multilocus.txt", 
              help="Multilocus risk weights.")
@click.option('-w',"--weights", type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.txt", help=OPTION_weights)
def sharp2019(geno, info, weights, multilocus):
	"""Calculate Gene Risk Score based on Sharp et al 2019. (DO NOT USE: UNDER DEVELOPMENT)."""
	check_args(geno, info, weights, multilocus)
	snpinfo = snps.ReadInfo(info)
	grs = riskscore.sharp2019(snpinfo, risks=weights, multirisks=multilocus)
	process_geno(geno, grs)

@main.command()
@click.option('--vcf', type=pkclick.gzFile(mode='rb'), help=OPTION_vcf)
@click.option('-m','--multilocus',  type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.multilocus.txt", 
              help="Multilocus risk weights.")
@click.option('-w','--weights', type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt", help=OPTION_weights)
def test(vcf, weights, multilocus):
	"""FOR TESTING PURPOSES ONLY; DO NOT USE!"""
	import vcf as vcf_
	vcfdata = snps.ReadVCF(vcf_.Reader(vcf), drop_genotypes=False)
	sbjgeno = transpose_geno([snp.genotype() for snp in vcfdata.values()])
	snpinfo = [snp.drop_genotype() for snp in vcfdata.values()]
	grs = riskscore.oram2016(snpinfo, risks=weights, multirisks=multilocus)
	for subjectid, gt in sbjgeno.items():
		print(subjectid + "\t" + str(grs.calc(gt)))

# --%%  END: Commands  %%--
#
##################################################




##################################################
#
# --%%  RUN: Subroutines  %%--

def check_args(geno, info, *args):
	if not (geno and info):
		sys.exit("ERROR: No input files specified, add both '--geno' and '--info' options.")
	return 1

# Another quick'n'dirty function which should probably be done cleaner
def transpose_geno(genoiter):
	import collections
	dict_dicts = collections.defaultdict(dict)
	for geno_dict in genoiter:
		for sbjid, gt in geno_dict.items():
			dict_dicts[sbjid][gt.ID] = gt
	return dict_dicts

# Should probably write an actual parser function which converts geno to instances of my genotype class...
def process_geno(geno, grs):
	for gtype in csv.DictReader(geno):
		subjectid = gtype.pop("")
		for k in gtype:
			gtype[k] = sum(float(x) for x in re.split("[/|]+", gtype[k])) 
		print(subjectid + "\t" + str(grs.calc(gtype)))


# --%%  END: Subroutines  %%--
#
##################################################


if __name__ == '__main__':
	main()


