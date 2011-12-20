# $Id$
#
# Copyright (C) 2001-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem
import collections

def _isCallable(thing):
  return (hasattr(collections,'Callable') and isinstance(thing,collections.Callable)) or \
              hasattr(thing,'__call__')   

_descList=[]
def _setupDescriptors(namespace):
  global _descList,descList
  from rdkit.Chem import GraphDescriptors,MolSurf,Lipinski,Fragments,Crippen
  from rdkit.Chem.EState import EState_VSA
  mods = [GraphDescriptors,MolSurf,EState_VSA,Lipinski,Crippen,Fragments]

  otherMods = [Chem]

  for nm,thing in namespace.iteritems():
    if nm[0]!='_' and _isCallable(thing):
      _descList.append((nm,thing))
  
  others = []
  for mod in otherMods:
    tmp = dir(mod)
    for name in tmp:
      if name[0] != '_':
        thing = getattr(mod,name)
        if _isCallable(thing):
          others.append(name)

  for mod in mods:
    tmp = dir(mod)

    for name in tmp:
      if name[0] != '_' and name[-1] != '_' and name not in others:
        # filter out python reference implementations:
        if name[:2]=='py' and name[2:] in tmp:
          continue
        thing = getattr(mod,name)
        if _isCallable(thing):
          namespace[name]=thing
          _descList.append((name,thing))
  descList=_descList

from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors
MolWt = lambda *x,**y:_rdMolDescriptors._CalcMolWt(*x,**y)
MolWt.version=_rdMolDescriptors._CalcMolWt_version
MolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> MolWt(Chem.MolFromSmiles('CC'))
  30.07...
  >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  53.49...

"""

HeavyAtomMolWt=lambda x:MolWt(x,True)
HeavyAtomMolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
  24.02...
  >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  49.46...

"""
HeavyAtomMolWt.version="1.0.0"

def NumValenceElectrons(mol):
  """ The number of valence electrons the molecule has

  >>> NumValenceElectrons(Chem.MolFromSmiles('CC'))
  14.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)O'))
  18.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)[O-]'))
  18.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)'))
  12.0
  

  """
  tbl = Chem.GetPeriodicTable()
  accum = 0.0
  for atom in mol.GetAtoms():
    accum += tbl.GetNOuterElecs(atom.GetAtomicNum())
    accum -= atom.GetFormalCharge()
    accum += atom.GetTotalNumHs()

  return accum
NumValenceElectrons.version="1.0.0"

def _oyMolecularFormula(mol):
   """ DEPRECATED: this is now implemented in C++
   Return the molecular formula

   contribution from Andrew Dalke
   
   """
   import collections

   # Count the different atom types
   counts = collections.defaultdict(int)
   charge = 0
   for atom in mol.GetAtoms():
       symb = atom.GetSymbol()
       counts[symb]+=1
       # Also track the implicit hydrogen counts and total charge
       num_hydrogens = atom.GetTotalNumHs()
       # A counts key should only exist if an element is present.
       # Do this test to prevent the creation of an "H": 0 for
       # things like [C], which have no implicit hydrogens.
       if num_hydrogens:
           counts['H'] += num_hydrogens
       charge += atom.GetFormalCharge()

   # Alphabetize elements by name
   elements = sorted(counts)

   # Put into Hill system order:
   # If there are carbons then they go first, followed by hydrogens,
   # then followed alphabetically by the other elements.
   # If there are no carbons then the elements are in alphabetical
   # order (including hydrogens)
   if "C" in elements:
       elements.remove("C")
       elements.insert(0, "C")
       if "H" in elements:
           elements.remove("H")
           elements.insert(1, "H")

   # Include the count, so {"C": 1, "H": 4} becomes ["C", "H4"]
   formula_terms = []
   for element in elements:
       formula_terms.append(element)
       if counts[element] > 1:
           formula_terms.append(str(counts[element]))

   # Handle the total charge. The result will be NH4+, Ca+2, etc.
   if charge == 0:
       pass
   elif charge == 1:
       formula_terms.append("+")
   elif charge == -1:
       formula_terms.append("-")
   else:
       formula_terms.append("%+d" % charge)

   return "".join(formula_terms)
MolecularFormula=lambda x:_rdMolDescriptors.CalcMolFormula(x)
MolecularFormula.version=_rdMolDescriptors._CalcMolFormula_version

_setupDescriptors(locals())



#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"],optionflags=doctest.ELLIPSIS)

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
