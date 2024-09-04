
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
recognized by having rsIDs separated by commas, while multiple genotypes can be separated with ':' or ';'. If the
interactions are more complicated, then try contacting the author (<fls530@ku.dk>) as it may be possible to implement a
specific command for the PGS. See the 'oram2016' and 'sharp2019' commands for examples of such.

Example of how to encode an interaction (from Oram2016):

\b
rsID                  effect_allele   effect_weight
rs2187668,rs7454108   C/T:T/C         3.87

""",
)

OPTION = namedtuple('Options', ['conflict','geno','info','log','logfile','n','pgs','vcf','weights','multiweights'])(
	conflict = """Handling of weight calculation when more than the maximum haplotypes are inferred. 'ignore' sets subject score to NA. 'Rank' uses the background frequencies from .. discarding the least common haplotypes.""",
	geno     = """Geno file of the type created by SNPextractor.""",
	info     = """Info file of the type created by SNPextractor.""",
	log      = """Control logging. Valid levels: 'debug', 'info', 'warning', 'error', 'critical'.""",
	logfile  = """Write log to FILE.""",
	n        = """The denominator to use in calculating the arithmetric mean of scores. Set to '1' to disable mean calculation.""",
	pgs      = """Risk score file obtained from the PGSCatalog. See: https://www.pgscatalog.org/.""",
	vcf      = """Read sample genotypes from VCF File as input.""",
	weights  = """Risk score file, possibly obtained from the PGSCatalog. See description below on 'Column-Based Datafiles'.\n""",
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
		self.epilog = EPILOG.fileformat + EPILOG.legal

class IntCommand(StdCommand):
	def __init__(self, *args, epilog=None, **kwargs):
		super().__init__(*args, **kwargs)
		self.epilog = EPILOG.fileformat + EPILOG.multiformat + EPILOG.legal


@click.group()
@click.version_option(version=__version__)
#@click.option('--dbsnp', help="what dbsnp should we use?")
@click.option('--log', default="warning", help=OPTION.log, show_default=True)
@click.option('--logfile', type=str, default=None, show_default="STDERR", help=OPTION.logfile)
def main(log, logfile):
	"""Calculate a Genetic Risk Score (GRS) for a list of subjects based on predefined weighted risks."""
	# Clear existing loggers...
	for handler in logging.root.handlers[:]:
		logging.root.removeHandler(handler)
	try:
		log_num = getattr(logging, log.upper())
	except AttributeError:
		raise ValueError(f"Invalid log level: '{log}'")
	logging.basicConfig(filename=logfile, filemode='w', level=log_num)
	logging.info(f"Logging set to level={log_num}, filename={logfile}")


@main.command(cls=IntCommand, no_args_is_help=True)
@click.option('-n','--denominator', type=click.FLOAT, help=OPTION.n)
@click.option('-w','--pgs','--weights', type=click.PGSFile(), required=True, help=OPTION.weights)
def aggregate(denominator, pgs, vcf):
	"""Calculate Aggregated (linear) Risk Score from user-provided weights."""
	from pkrs import RiskScore
	rs = RiskScore.FromPGS(pgs, N=denominator)
	calc_and_report(rs, vcf)


oram2016_weights_default = os.environ.get('RISKSCORE_ORAM2016_WEIGHTS', None)
@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-n','--denominator', type=click.FLOAT, default=58, show_default=True, help=OPTION.n)
@click.option('-w', '--pgs', '--weights', type=click.PGSFile(), envvar='RISKSCORE_ORAM2016_WEIGHTS', default=oram2016_weights_default, required=True, help=OPTION.weights, show_default=True)
def oram2016(denominator, pgs, vcf):
	"""Calculate T1D Risk Score based on Oram et al 2016.

\b
A type 1 diabetes genetic risk score can aid discrimination between type 1 and
type 2 diabetes in young adults.
RA Oram, K Patel, A Hill, B Shields, TJ McDonald, A Jones, AT Hattersley,
MN Weedon.
Diabetes care 39 (3), 337-344.
https://doi.org/10.2337/dc15-1111
"""
	from pkrs import Oram2016
	rs = Oram2016.FromPGS(pgs, N=denominator, interaction_func=max)
	calc_and_report(rs, vcf)


sharp2019_weights_default = os.environ.get('RISKSCORE_SHARP2019_WEIGHTS', None)
@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--conflict', type=click.Choice(['Max','Mean','Ignore','Rank'], case_sensitive=False), default="Ignore", show_default=True, help=OPTION.conflict)
@click.option('-n', '--denominator', type=click.FLOAT, default=1, show_default=True, help=OPTION.n)
@click.option('-w', '--pgs', '--weights', type=click.PGSFile(), envvar='RISKSCORE_SHARP2019_WEIGHTS', default=sharp2019_weights_default, required=True, help=OPTION.weights, show_default=True)
def sharp2019(conflict, denominator, pgs, vcf):
	"""Calculate T1D Risk Score based on Sharp et al 2019.

\b
Development and standardization of an improved type 1 diabetes genetic risk
score for use in newborn screening and incident diagnosis.
SA Sharp, SS Rich, AR Wood, SE Jones, RN Beaumont, JW Harrison, DA Schneider,
JM Locke, JT, MN Weedon, WA Hagopian, RA Oram.
Diabetes Care 2019 Feb; 42(2): 200-207.
https://doi.org/10.2337/dc18-1785
"""
	logging.debug(f" Haplotype conflict resolution = {conflict}")
	if (conflict == "Max"):
		int_func = max
		hap_func = lambda x: sum(sorted(x, reverse=True)[0:2])
	elif (conflict == "Ignore"):
		int_func = lambda x: "NA" if len(x) > 1 else sum(x)
		hap_func = lambda x: "NA" if len(x) > 2 else sum(x)
	elif (conflict == "Mean"):
		int_func = lambda x: sum(x) / len(x)
		hap_func = lambda x: sum(x) if len(x) <= 1 else 2 * sum(x) / len(x)
	elif (conflict == "Rank"):
		# Setting up the sorting - The index from table below is what we should sort on (Though we will need to round the shit to avoid real mis/match snafus)
		klitz_sort = lambda x: [2.31, -0.24, 3.63, -1.94, 0.05, 0.17, -0.51, 0.17, 2.16, -0.69, 1.08, -0.65, 1.06, -1.47, 3.14, -0.87, 0.61, -1.15].index(round(x, 2))
		int_func = lambda x: sorted(x, key=klitz_sort)[0]
		hap_func = lambda x: sum(sorted(x, reverse=True)[0:2]) # Yes, by chance the frequencies and beta weights correlate, so highest is most frequent

	from pkrs import Sharp2019
	rs = Sharp2019.FromPGS(pgs, N=denominator, interaction_func=int_func, haplotype_func=hap_func)
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
		logging.debug(f" Processing sample={sample} with alleles={alleles}")
		score = rs.calc(alleles)
		logging.debug(f" Sample = {sample}; Total Score = {score}")
		print(f"{sample}\t{score}")

# --%%  END: Subroutines  %%--
#
##################################################

