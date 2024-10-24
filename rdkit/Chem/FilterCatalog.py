#  Copyright (C) 2015 Novartis Institute for BioMedical Research
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import sys

from rdkit import Chem, rdBase
from rdkit.Chem.rdfiltercatalog import *
MatchTypeVect = rdBase.MatchTypeVect

class FilterMatcher(PythonFilterMatcher):
  """FilterMatcher - This class allows creation of Python based
    filters.  Subclass this class to create a Filter useable
    in a FilterCatalogEntry

    Simple Example:

    from rdkit.Chem import rdMolDescriptors
    class MWFilter(FilterMatcher):
      def __init__(self, minMw, maxMw):
          FilterMatcher.__init__(self, "MW violation")
          self.minMw = minMw
          self.maxMw = maxMw

      def IsValid(self):
         return True

      def HasMatch(self, mol):
         mw = rdMolDescriptors.CalcExactMolWt(mol)
         return not self.minMw <= mw <= self.maxMw
    """

  def __init__(self, name="Unamed FilterMatcher"):
    self.name = name
    PythonFilterMatcher.__init__(self, self)

  def HasMatch(self, mol):
    """Return True if the filter matches the molecule"""
    raise NotImplementedError("Need to implement HasMatch(mol) in a subclass of %s",
                              self.__class__.__name__)

  def GetMatches(self, mol, matchVect):
    """Return True if the filter matches the molecule
        (By default, this calls HasMatch and does not modify matchVect)
        
        matchVect is a vector of FilterMatch's which hold the matching
        filter and the matched query_atom, mol_atom pairs if applicable.
        To append to this vector:
        v = MatchTypeVect()
        v.append(IntPair( query_atom_idx, mol_atom_idx ) )
        match = FilterMatch(self, v)
        matchVect.append( match )
        """
    return self.HasMatch(mol)

  def IsValid(self, mol):
    """Must override this function"""
    raise NotImplementedError("IsValid must be implemented in a subclass of %s",
                              self.__class__.__name__)

  def GetName(self):
    return self.name
