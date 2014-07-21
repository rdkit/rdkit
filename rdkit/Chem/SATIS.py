# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""  Functionality for SATIS typing atoms

"""
from __future__ import print_function
from rdkit import Chem
from rdkit.six.moves import xrange

_debug = 0

#
#  These are SMARTS patterns for the special cases used in
#  SATIS typing.
#
aldehydePatt = Chem.MolFromSmarts('[CD2]=[OD1]')
ketonePatt = Chem.MolFromSmarts('[CD3]=[OD1]')
amidePatt = Chem.MolFromSmarts('[CD3](=[OD1])-[#7]')
esterPatt = Chem.MolFromSmarts('C(=[OD1])-O-[#6]')
carboxylatePatt = Chem.MolFromSmarts('C(=[OD1])-[OX1]')
carboxylPatt = Chem.MolFromSmarts('C(=[OD1])-[OX2]')

specialCases = ((carboxylatePatt,97),
                (esterPatt,96),
                (carboxylPatt,98),
                (amidePatt,95),
                (ketonePatt,94),
                (aldehydePatt,93))


def SATISTypes(mol,neighborsToInclude=4):
  """ returns SATIS codes for all atoms in a molecule

   The SATIS definition used is from:
   J. Chem. Inf. Comput. Sci. _39_ 751-757 (1999)

   each SATIS code is a string consisting of _neighborsToInclude_ + 1
   2 digit numbers

   **Arguments**

     - mol: a molecule

     - neighborsToInclude (optional): the number of neighbors to include
       in the SATIS codes

   **Returns**

     a list of strings nAtoms long

  """
  global specialCases
  
  nAtoms = mol.GetNumAtoms()
  atomicNums = [0]*nAtoms
  atoms = mol.GetAtoms()
  for i in xrange(nAtoms):
    atomicNums[i] = atoms[i].GetAtomicNum()

  nSpecialCases = len(specialCases)
  specialCaseMatches = [None]*nSpecialCases
  for i,(patt,idx) in enumerate(specialCases):
    if mol.HasSubstructMatch(patt):
      specialCaseMatches[i] = mol.GetSubstructMatches(patt)
    else:
      specialCaseMatches[i] = ()
    
  codes = [None]*nAtoms
  for i in range(nAtoms):
    code = [99]*(neighborsToInclude+1)
    atom = atoms[i]
    atomIdx = atom.GetIdx()
    code[0] = min(atom.GetAtomicNum(),99)
    bonds = atom.GetBonds()
    nBonds = len(bonds)
    otherIndices = [-1]*nBonds
    if _debug: print(code[0],end='')
    for j in range(nBonds):
      otherIndices[j] = bonds[j].GetOtherAtom(atom).GetIdx()
      if _debug: print(otherIndices[j],end='')
    if _debug: print()  
    otherNums = [atomicNums[x] for x in otherIndices] + \
                [1]*atom.GetTotalNumHs()
    otherNums.sort()

    nOthers = len(otherNums)
    if nOthers > neighborsToInclude:
      otherNums.reverse()
      otherNums = otherNums[:neighborsToInclude]
      otherNums.reverse()
      for j in range(neighborsToInclude):
        code[j+1] = min(otherNums[j],99)
    else:
      for j in range(nOthers):
        code[j+1] = min(otherNums[j],99)
      if nOthers < neighborsToInclude and code[0] in [6,8]:
        found = 0
        for j in range(nSpecialCases):
          for matchTuple in specialCaseMatches[j]:
            if atomIdx in matchTuple:
              code[-1] = specialCases[j][1]
              found = 1
              break
          if found:
            break
        
    codes[i] = ''.join(['%02d'%(x) for x in code])
  return codes

if __name__ == '__main__':
  smis = ['CC(=O)NC','CP(F)(Cl)(Br)(O)',
          'O=CC(=O)C','C(=O)OCC(=O)O','C(=O)[O-]']
  for smi in smis:
    print(smi)
    m = Chem.MolFromSmiles(smi)
    codes = SATISTypes(m)
    print(codes)
