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

Version = "0.1 (Development Version)"

import argparse
import pkcsv as csv
import pksnps as snps
import pkrs as riskscore
import sys
 
# --%%  END: Perform Basic Setup  %%--
#
##################################################



##################################################
#
# --%%  RUN: Subroutines  %%--

#def showversion():
#	print(Version)
#	sys.exit(0)

def get_algorithm(args):
	snpinfo = snps.ReadInfo(args.info)
	if args.oram2016:
		grs = riskscore.oram2016(snpinfo)
	else:
		sys.exit("Algorithm Error!")
	return grs



##################################################
#
# --%%  RUN: Main Program  %%--

# Setup ArgParse Here...

def setup_args():
	parser = argparse.ArgumentParser(description="THIS SCRIPT IS STILL EXPERIMENTAL; USE WITH CAUTION",
	                         epilog="Fisk")
	parser.add_argument('-g', '--geno', default=None,
		help="Geno file")
	parser.add_argument('-i', '--info', default=None,
		help="Info file")
	parser.add_argument('--oram2016', const="risk_score/oram2016.weights_hla.txt", default="", nargs="?",
		help="Oram et al 2016")
	parser.add_argument('-v', '--version', action='version', version=Version, help="Print Version and exit.")
# choice of algorithm
# optional weights file
	args = parser.parse_args()
	if not (args.geno or args.info):
		parser.error('No input files specified, add --geno or --info (or both)')
	args.info = args.info if args.info else [args.geno.replace(".geno.",".info.")]
#	if args.version:
#		showversion()
	return args


# Execute Main Program from Here

def main(args):
	grs = get_algorithm(args)

# loop over input
	with open(args.geno) as f:
		for geno in csv.DictReader(f):
			subjectid = geno.pop("")
			print(subjectid + "\t" + str(grs.calc(geno)))
	return 0

if __name__ == '__main__':
	args = setup_args()
	state = main(args)
	sys.exit(state)

# --%%  END: Main program  %%--
#
##################################################





