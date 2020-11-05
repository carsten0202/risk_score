#!/home/fls530/anaconda3/bin/python
 
import click
import gzip
import io

magic_dict = {
	 b"\x1f\x8b\x08": "gz",
	 b"\x42\x5a\x68": "bz2",
	 b"\x50\x4b\x03\x04": "zip"
}

class gzFile(click.File):
	def convert(self, value, param, ctx):
		f = super().convert(value, param, ctx)
		if self._getziptype(f) is not None:
			return io.TextIOWrapper(gzip.GzipFile(fileobj=f))
		return f

	@staticmethod
	def _getziptype(f):
		if f.seekable():
			file_start = f.read(max(len(x) for x in magic_dict))
			f.seek(0)
			for magic, filetype in magic_dict.items():
				if file_start.startswith(magic):
					return filetype
		return None

