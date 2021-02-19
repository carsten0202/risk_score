#!/home/fls530/anaconda3/bin/python
 
import click

#
# -%  Class pkclick.gzFile  %-

class gzFile(click.File):
	"""A Class for detecting compressed files and automagically decrompress them."""
	magic_dict = {
		b"\x1f\x8b\x08"     : "gz",
		b"\x42\x5a\x68"     : "bz2",
		b"\x50\x4b\x03\x04" : "zip"
	}

	def convert(self, value, param, ctx):
		import gzip
		import io
		f = super().convert(value, param, ctx)
		try:
			if self._getziptype(f, self.magic_dict) is not None:
				return io.TextIOWrapper(gzip.GzipFile(fileobj=f))
		except UnicodeDecodeError:
			self.fail("Could not interpret input. Did you remember to use binary mode? eg gzFile(mode='rb')")
		return f

	@staticmethod
	def _getziptype(f, magic_dict):
		if f.seekable():
			file_start = f.read(max(len(x) for x in magic_dict))
			f.seek(0)
			for magic, filetype in magic_dict.items():
				if file_start.startswith(magic):
					return filetype
		return None



#
# -%  CLASS: pkclick.csv  %-

class CSV(click.ParamType):
	"""A Class for handling option values that are comma separated values."""
	name = "csv_list"
	def convert(self, value, param, ctx):
		import csv
		value = super().convert(value, param, ctx)
		return next(csv.reader([value]))



#
# -%  CLASS: Samples  %-

class SampleList(gzFile):
	"""Obtain a list of samples from a file (or '-'... maybe?)."""
	def convert(self, value, param, ctx):
		f = super().convert(value, param, ctx)
		if self._isVCF(f):
			try:
				import vcf
				vcf_r = vcf.Reader(f, compressed=False)
				record = next(vcf_r)
				return [subject.sample for subject in record.samples]
			except ImportError:
				self.fail("Encountered VCF imput but could not find the PyVCF module. Either install with 'pip install PyVCF' or provide different input.")
		self.fail("Could not parse input. Please try a different file format.")

	@staticmethod
	def _isVCF(f):
		if f.seekable():
			file_start = f.readline()
			f.seek(0)
			return file_start.startswith("##fileformat=VCFv4")
		return None

