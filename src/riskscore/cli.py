
##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

# import click
from collections import namedtuple
import logging
import os
import sys

import pklib.pkclick as click
from riskscore.version import __version__

EPILOG = namedtuple('Options', ['fileformat','multiformat','legal'])(
fileformat = """

VARIANT MATCHING:

Variants are matched based on rsIDs and nucleotides. PGS files and vcf files must be formatted with rsIDs or proper
matching of variants will not occur.

COLUMN-BASED DATAFILES:

Accepted input is either a PGS accession, which is then downloaded from the catalog, or a tabular file formatted
according to the PGS standard. Column separators are auto-detected from the input and should work for tab-separated
files, comma-seperated (csv) files, and space-separated files. Columns with unrecognized names are ignored.

\b
https://www.pgscatalog.org/downloads/#dl_scoring_files

""",
legal = """

LEGAL:

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Written by Carsten Friis Rundsten <fls530@ku.dk>
""",
multiformat = """

INTERACTIONS:

The PGS catalog file format is inconsistent when describing weighted interactions. Limited support is provided here if
weights are appropriately formatted, which may require the user to modify the PGS file after download. Interactions are
recognized by having rsIDs separated by commas, while multiple genotypes can be separated with ':' or ';'.

Example of an interaction from Oram2016:

\b
rsID                  effect_allele   effect_weight
rs2187668,rs7454108   C/T:T/C         3.87

""",
)

OPTION = namedtuple('Options', ['geno','info','log','n','pgs','vcf','weights','multiweights'])(
	geno = """Geno file of the type created by SNPextractor.""",
	info = """Info file of the type created by SNPextractor.""",
	log  = """Control logging. Valid levels: 'debug', 'info', 'warning', 'error', 'critical'.""",
	n    = """The denominator to use in calculating the arithmetric mean of scores. Set to '1' to disable mean calculation.""",
	pgs  = """Risk score file obtained from the PGSCatalog. See: https://www.pgscatalog.org/.""",
	vcf  = """Read sample genotypes from VCF File as input.""",
	weights = """Risk score file, possibly obtained from the PGSCatalog. See description below on 'Column-Based Datafiles'.\n""",
	multiweights = """Multi-locus risk weights file. See format description below on 'Column-Based Datafiles'.\n"""
)

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  RUN: Define Commands  %%--

class StdCommand(click.Command):
	def __init__(self, *args, epilog=None, **kwargs):
		super().__init__(*args, **kwargs)
		self.params.insert(0, click.Option(('--vcf',), type=click.VCFFile(), help=OPTION.vcf))
#		self.params.insert(0, click.Option(('-i','--info',), type=click.File(), help=OPTION.info))
#		self.params.insert(0, click.Option(('-g','--geno',), type=click.File(), help=OPTION.geno))
		self.epilog = EPILOG.fileformat + EPILOG.multiformat + EPILOG.legal

@click.group()
@click.version_option(version=__version__)
@click.option('--dbsnp', help="what dbsnp should we use?")
@click.option('--log', default="warning", help=OPTION.log, show_default=True)
def main(dbsnp, log):
	"""Calculate a Genetic Risk Score (GRS) for a list of subjects based on predefined weighted risks."""
	# Clear existing loggers...
	for handler in logging.root.handlers[:]:
		logging.root.removeHandler(handler)
	try:
		log_num = getattr(logging, log.upper())
	except AttributeError:
		raise ValueError(f"Invalid log level: '{log}'")
	logging.basicConfig(level=log_num)


@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-n','--denominator', type=click.FLOAT, help=OPTION.n)
@click.option('-w','--pgs','--weights', type=click.PGSFile(), default=None, help=OPTION.weights)
def aggregate(denominator, pgs, vcf):
	"""Calculate Aggregated (linear) Risk Score from user-provided weights."""
	assert pgs, f"ERROR: You must specify a weights file."
	from pkrs import RiskScore
	rs = RiskScore.FromPGS(pgs, N=denominator)
	calc_and_report(rs, vcf)


oram2016_weights_default = os.environ.get('RISKSCORE_ORAM2016_WEIGHTS', None)
@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-n','--denominator', type=click.FLOAT, default=58, show_default=True, help=OPTION.n)
@click.option('-w', '--pgs', '--weights', type=click.PGSFile(), envvar='RISKSCORE_ORAM2016_WEIGHTS', default=oram2016_weights_default, required=True, help=OPTION.weights, show_default=True)
def oram2016(denominator, pgs, vcf):
	"""Calculate Gene Risk Score based on Oram et al 2016.

\b
A type 1 diabetes genetic risk score can aid discrimination between type 1 and
type 2 diabetes in young adults.
RA Oram, K Patel, A Hill, B Shields, TJ McDonald, A Jones, AT Hattersley,
MN Weedon.
Diabetes care 39 (3), 337-344.
https://doi.org/10.2337/dc15-1111
"""
	from pkrs import Interaction as Oram2016
	rs = Oram2016.FromPGS(pgs, N=denominator, interaction_func=max)
	calc_and_report(rs, vcf)


sharp2019_weights_default = os.environ.get('RISKSCORE_SHARP2019_WEIGHTS', None)
@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-n','--denominator', type=click.FLOAT, default=1, show_default=True, help=OPTION.n)
@click.option('-w', '--pgs', '--weights', type=click.PGSFile(), envvar='RISKSCORE_SHARP2019_WEIGHTS', default=sharp2019_weights_default, required=True, help=OPTION.weights, show_default=True)
def sharp2019(denominator, pgs, vcf):
	"""Calculate Gene Risk Score based on Sharp et al 2019.

\b
Development and standardization of an improved type 1 diabetes genetic risk
score for use in newborn screening and incident diagnosis.
SA Sharp, SS Rich, AR Wood, SE Jones, RN Beaumont, JW Harrison, DA Schneider,
JM Locke, JT, MN Weedon, WA Hagopian, RA Oram.
Diabetes Care 2019 Feb; 42(2): 200-207.
https://doi.org/10.2337/dc18-1785
"""
	from pkrs import Sharp2019
	rs = Sharp2019.FromPGS(pgs, N=denominator, interaction_func=max)
	calc_and_report(rs, vcf)


@main.command(no_args_is_help=True, hidden=True)
@click.option('--vcf', type=click.File(), default="/home/fls530/dev/testdata/oram_sharp.sorted.vcf.gz", help=OPTION.vcf)
@click.option('-w', '--weights', '--pgs', type=str, default="/home/fls530/dev/riskscore/weights/oram2016.txt.gz", help=OPTION.weights, show_default=True)
def test(vcf, weights):
	"""FOR TESTING PURPOSES ONLY; DO NOT USE!"""
	logging.warning(f"FOR DEVELOPMENTAL TESTING PURPOSES ONLY; DO NOT USE!")

	from pkrs import Interaction
	from pksnp import PopulationAlleles

	import pgscatalog.core
	pgs = pgscatalog.core.ScoringFile(weights)
	try:
		pgs.download(".") # type: ignore
	except FileExistsError:
		logging.debug(f": File exists - Download aborted")

	grs = Interaction.FromPGS(pgs)
	population = PopulationAlleles(vcf)

	for sample, alleles in population.items():
		logging.info(f" Processing sample={sample}")
		logging.debug(f" ...with alleles={alleles}")
		score = grs.calc(alleles)
		logging.debug(f" Sample = {sample}; Total Score = {score}")
		print(f"{sample}\t{score}")

# --%%  END: Define Commands  %%--
#
##################################################




##################################################
#
# --%%  RUN: Subroutines  %%--

def check_args(geno=None, info=None, vcf=None, *args):
	"""Check that EITHER Geno & Info pair or VCF file is given."""
	if bool(vcf) ^ bool(geno and info):
		return vcf or geno
	sys.exit("ERROR: No correct input files specified, use EITHER '--vcf' OR BOTH '--geno' and '--info' options. Use '--help' for help.")

def calc_and_report(rs, vcf):
	"""The standardised part of calculating"""
	from pksnp import PopulationAlleles
	population = PopulationAlleles(vcf, filter_ids=rs.rsids)
	for sample, alleles in population.items():
		logging.info(f" Processing sample={sample}")
		logging.debug(f" ...with alleles={alleles}")
		score = rs.calc(alleles)
		logging.debug(f" Sample = {sample}; Total Score = {score}")
		print(f"{sample}\t{score}")



# --%%  END: Subroutines  %%--
#
##################################################

