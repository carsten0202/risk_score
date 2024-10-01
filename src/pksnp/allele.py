###########################################################
#
# ---%%%  allele.py: Module to hold hashable alleles for matching PRS scores  %%%---
#

import logging

logger = logging.getLogger(__name__)


# This is intended as a hashable allele to be used as key in a {allele:call} dict structure
        
class Allele:
    """
    Class to hold an Allele. Can be used as key in e.g. dict(key) = (call) or dict(key) = risk_weight
	Must be immutable; DO NOT change after __init__

    chromosome (str): The chromosome
    position (int): Chromosomal position
    rsid (str): The rsid of the variant
    allele (str): The allele (nucleotides)
    build (str): Genomic build (see _validate)
    """
    def __init__(self, chromosome=None, position=None, rsid=None, allele=None, build=None):
        assert allele,  f"Allele cannot be empty or NoneType."
        self.chromosome = str(chromosome)
        self.position = int(position) if position is not None else None
        self.rsid = str(rsid)
        self.allele = str(allele)
        self.build = build
        Allele._validate(self)

    def __eq__(self, other):
        # return (self.chromosome == other.chromosome and
        #         self.position == other.position and
        #         self.rsid == other.rsid and
        #         self.allele == other.allele)
        if isinstance(other, Allele):
            return (self.rsid == other.rsid and
                    self.allele == other.allele)
        return NotImplemented

    def __hash__(self):
#        return hash((self.chromosome, self.position, self.rsid, self.allele))
        return hash((self.rsid, self.allele))

    def __repr__(self):
        return f"Allele({self.rsid}, {self.allele})"

    @staticmethod
    def _validate(obj):
        builds = [None, "GRCh37", "GRCh38"]
        assert obj.build in builds, f"Build type ({obj.build}) must be one of {builds}."
        assert obj.rsid or (obj.chromosome and obj.position), "Allele must have either a chr:pos or an rsid."

    @property
    def chrom(self):
        """Abbrevation of chromosome."""
        return self.chromosome

    @property
    def pos(self):
        """Abbrevation of position."""
        return self.position
