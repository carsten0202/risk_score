
__version__ = """2.1.0"""
# v1.5.0: Fixed Sharp2019. Kinda. Its still wonky, and has trouble matching the alleles properly. But it works with some hand-holding.
# v1.6.0: Fixed even more bugs. Sharp is a little less wonky. Still has trouble matching the alleles.
# v1.6.1: Bugfix
# v1.6.2: Integrated the nested_lookup functions in sharp and multirisk into one generic function

# v2.0.0: We'll make v2 using pgscatalog tools, so, yes, it's v2. (Even if the tools works like shit...)
# v2.1.0: Added option to control haplotype conflicts in Sharp2019

# Notes and TODOs:
# TODO: Could integrate the database code from snptool to speedup translation between rsid and pos
# TODO: When reading from VCF we should include an option for selecting the field you want; eg GT or DS, etc.
# TODO: Where should we put score files downloaded with pgscatalog tools? Currently they end in '.', which may be ok?
