
__version__ = """1.6.2"""

# v1.5.0: Fixed Sharp2019. Kinda. Its still wonky, and has trouble matching the alleles properly. But it works with some hand-holding.
# v1.6.0: Fixed even more bugs. Sharp is a little less wonky. Still has trouble matching the alleles.
# v1.6.1: Bugfix
# v1.6.2: Integrated the nested_lookup functions in sharp and multirisk into one generic function

# Notes and TODOs:

# TODO: When reading from VCF we should include an option for selecting the field you want; eg GT or DS, etc.
# TODO: The filtering isn't complete in pkcsv. I only filter '#' at the beginning, but we should do it in the whole file.

