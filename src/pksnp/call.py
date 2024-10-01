###########################################################
#
# ---%%%  call.py: Module to hold variant calls  %%%---
#

import logging

logger = logging.getLogger(__name__)


# There is no proper support for dosage scores in the script yet, but we could add it. Will likely need a dosage class:
        
class Call(object):
    """
    Class to hold a simple Variant Call. Introduced to get consistent behaviour regardless of whether calls are
    genotypes, dosages or probabilities. Can be used as value in SampleAlleles dicts e.g. dict(Allele:call).

    Attributes:
    _store (list): List of ploidy length with the probabilities for the allele being present in each copy.
    """
    def __init__(self, index:int, DS=None, GT=None, GP=None, HDS=None):
        """
        index (int):       Index of allele under consideration. This is 0 for reference, 1 for first alt, 2 for second, etc.
        DS (float):        Dosage score
        GT (list(int)):    List of integer genotype calls
        GP (NotDone!):     Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1
        HDS (list(float)): List of  Haploid Alternate Allele Dosage scores
        """
        if HDS:
            assert index in [0,1], "Spurious dosage scores for multivariant SNP encountered."
            self._store = HDS if index else [abs(ds-1) for ds in HDS]
        elif GT:
            self._store = [int(x==index) for x in GT]
        elif GP:
            assert False, "GP - Not supported yet!"
            self._store = []
        elif DS:
            assert index in [0,1], "Spurious dosage scores for multivariant SNP encountered."
            self._store = [min(1,DS), max(0,DS-1)] if index else [min(1,abs(DS-2)), max(0,abs(DS-2)-1)]

    def __add__(self, value, /):
        return self.dosage + value

    def __bool__(self, /):
        return any([x>0.5 for x in self.hds])

    def __eq__(self, other):
        if isinstance(other, Call):
            return self.genotype == other.genotype
        else:
            return self.genotype == other

    def __ge__(self, value, /):
        return sum(self.genotype) >= value

    def __gt__(self, value, /):
        return sum(self.genotype) > value

    def __le__(self, value, /):
        return not self > value

    def __lt__(self, value, /):
        return not self >= value

    def __ne__(self, value, /):
        return not self == value

    def __iter__(self):
        yield self.dosage

    __radd__ = __add__

    def __repr__(self):
        return f"Call({self.hds})"


    @property
    def dosage(self):
        """The Call expressed as a dosage score. Note that this assumes the call is the alternate as the Call class has
        no knowledge of true refs and alts."""
        return sum(self.hds)

    @property
    def genotype(self):
        """The Call expressed as list of maximum likelihood genotypes"""
        return [int(x>0.5) for x in self.hds]

    @property
    def hds(self):
        """The Call expressed as Haploid Allele Dosage. Note that this assumes the call is the alternate as the Call
        class has no knowledge of true refs and alts."""
        return self._store

    @property
    def is_heterogeneous(self):
        return bool(self.genotype[0]) != bool(self.genotype[1])

    @property
    def is_homogeneous(self):
        return all([x>0.5 for x in self.hds])

# These are not finished
# __sub__(self, value, /)
#     Return self-value.

# __rsub__(self, value, /)
#     Return value-self.

# __hash__(self, /)
#     Return hash(self).

