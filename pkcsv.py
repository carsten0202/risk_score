#!/home/fls530/anaconda3/envs/snakemake/bin/python

from csv import DictWriter, writer
import csv

def reader(f, dialect=None, **fmtparams):
	line = f.readline()
	dialect = csv.Sniffer().sniff(line)
	rows = next(csv.reader([line], dialect))
	rfunc = csv.reader(f, dialect=dialect, **fmtparams)
	return rfunc, rows

def DictReader(f, fieldnames=None, dialect=None, *args, **kwds):
	line = f.readline()
	dialect = csv.Sniffer().sniff(line)
	rows = next(csv.reader([line], dialect))
	rfunc = csv.DictReader(f, fieldnames=rows, dialect=dialect, *args, **kwds)
	return rfunc


