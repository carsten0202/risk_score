#!/home/fls530/anaconda3/bin/python

#
#  risk score ???
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



##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

Version = "0.4 (Development Version)"

import click
import pathlib
import pkcsv as csv
import pksnps as snps
import pkrs as riskscore
import re
import sys
 
# --%%  END: Perform Basic Setup  %%--
#
##################################################


EPILOG_FileFormat = """
Column-Based Datafiles:
Several options accept files containing data in tables/columns. Columns separators are auto-detected from the input and should work for tab-separated files, comma-seperated (csv) files, and space separated files. 
The file must have one and only one headline as columns are identified by name. The column order is ignored.

Any column name not recognized is ignored. Recognized names are:
ALLELE    - Identification of the weighted allele. Usually given as a nucleotide. 
BETA      - The weight to be used. Used unmodified. Takes precedent over 'ODDSRATIO'.
CHROM     - Chromosome name or number.
ODDSRATIO - If given without a 'BETA' column, then the natural logarithm of these values will be used as weights. If 'BETA' is given, then this column is ignored.
POS       - Chromosomal Position.
POSID     - An identifier of the type CHR:POS.

"""


##################################################
#
# --%%  RUN: Commands  %%--

class StdCommand(click.Command):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.params.insert(0, click.Option(('-i','--info',), type=click.File(), help='Info file.'))
		self.params.insert(0, click.Option(('-g','--geno',), type=click.File(), help='Geno file.'))
		self.epilog = "Test:" + EPILOG_FileFormat

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
	snpinfo = snp.ReadInfo(info)
	grs = riskscore.RiskScore(snpinfo, risks=weights)
	process_geno(geno, grs)

@main.command(cls=StdCommand)
@click.option('-m','--multilocus',  type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.multilocus.txt", 
              help="Multilocus risk weights.")
@click.option('-w','--weights', type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/oram2016.weights.txt", 
              help="Single locus risk weights file.")
def oram2016(geno, info, weights, multilocus):
	"""Calculate Gene Risk Score based on Oram et al 2016."""
	check_args(geno, info, weights, multilocus)
	snpinfo = snps.ReadInfo(info)
	grs = riskscore.oram2016(snpinfo, risks=weights, multirisks=multilocus)
	process_geno(geno, grs)

@main.command(cls=StdCommand)
@click.option('-m','--multilocus',  type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.multilocus.txt", 
              help="Multilocus risk weights.")
@click.option('-w',"--weights", type=click.File(), default=str(pathlib.Path(__file__).resolve().parent.absolute()) + "/sharp2019.weights.txt", 
              help='Single locus risk weights file.')
def sharp2019(geno, info, weights, multilocus):
	"""Calculate Gene Risk Score based on Sharp et al 2019. (DO NOT USE: UNDER DEVELOPMENT)."""
	check_args(geno, info, weights, multilocus)
	snpinfo = snps.ReadInfo(info)
	grs = riskscore.sharp2019(snpinfo, risks=weights, multirisks=multilocus)
	process_geno(geno, grs)

@main.command()
@click.option('--vcf', type=click.File(), help="Load VCF File")
def test(vcf):
	"""For testing purposes; Do Not Use!"""
	import vcf as vcf_
	vcfiter = vcf_.Reader(vcf)

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


