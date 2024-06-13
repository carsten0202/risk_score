
__version__ = """2.0.0"""
# v1.5.0: Fixed Sharp2019. Kinda. Its still wonky, and has trouble matching the alleles properly. But it works with some hand-holding.
# v1.6.0: Fixed even more bugs. Sharp is a little less wonky. Still has trouble matching the alleles.
# v1.6.1: Bugfix
# v1.6.2: Integrated the nested_lookup functions in sharp and multirisk into one generic function

# v1.7? v2.0? We will make some big changes for sharp2019 although the rest will be largely untouched.
#   We could fix some of the loading problems, and make the snp matching better...
# v2.0.0: We'll make v2 using pygscatalog, so, yes, it's v2.
#   And no: pygscatalog still works like shit! We'll stick with v2, but do it with help from chatGPT
#   We will need two calc functions. Yes, two.
#       One: The classic 'calc by sample' done in the grs classes. Slow as shit but can handle interactions
#       Two: Novel by snp with increments for each sample. Much faster but linear weights only.
#            This must be done in a VCF object or similar, since it's done while reading the VCF.
#       We do not need to implement two just yet.

# Notes and TODOs:
# TODO: Should integrate the database code from snptool to speedup translation between rsid and pos
#   TODO: The dbsnp thing is hardlined in right now... see 'sample_alleles'.
# TODO: Filtering PopAlleles so we only load the weighted ones?
#	Step one would require the grs to return all weighted rsids?
# TODO: When reading from VCF we should include an option for selecting the field you want; eg GT or DS, etc.
# TODO: The filtering isn't complete in pkcsv. I only filter '#' at the beginning, but we should do it in the whole file.
# TODO: Where should we put score files downloaded with pgstools?
