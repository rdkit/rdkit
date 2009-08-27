## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2006-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Import all RDKit chemistry modules
 
"""
from rdkit import rdBase
from rdkit import RDConfig
from rdkit.Geometry import rdGeometry
from rdkit.Chem import *
from rdDepictor import *
from rdChemReactions import *
import numpy

def EnumerateLibraryFromReaction(reaction,sidechainSets) :
  """ Returns a generator for the virtual library defined by
   a reaction and a sequence of sidechain sets

  >>> from rdkit import Chem
  >>> from rdkit.Chem import AllChem
  >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
  >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
  >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
  >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
  ['CNC=O', 'CCNC=O', 'CNC(=O)C', 'CCNC(=O)C']

  Note that this is all done in a lazy manner, so "infinitely" large libraries can
  be done without worrying about running out of memory. Your patience will run out first:

  Define a set of 10000 amines:
  >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))

  ... a set of 10000 acids
  >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))

  ... now the virtual library (1e8 compounds in principle):
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])

  ... look at the first 4 compounds:
  >>> [Chem.MolToSmiles(r.next()[0]) for x in range(4)]
  ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']

  
  """
  if len(sidechainSets) != reaction.GetNumReactantTemplates():
    raise ValueError,'%d sidechains provided, %d required'%(len(sidechainSets),reaction.getNumReactantTemplates())

  def _combiEnumerator(items,depth=0):
    for item in items[depth]:
      if depth+1 < len(items):
        v = _combiEnumerator(items,depth+1)
        for entry in v:
          l=[item]
          l.extend(entry)
          yield l
      else:
        yield [item]
  for chains in _combiEnumerator(sidechainSets):
    prodSets = reaction.RunReactants(chains)
    for prods in prodSets:
      yield prods

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
