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

import itertools

from rdkit import Chem

#  These are SMARTS patterns for the special cases used in SATIS typing.
aldehydePatt = Chem.MolFromSmarts('[CD2]=[OD1]')
ketonePatt = Chem.MolFromSmarts('[CD3]=[OD1]')
amidePatt = Chem.MolFromSmarts('[CD3](=[OD1])-[#7]')
esterPatt = Chem.MolFromSmarts('C(=[OD1])-O-[#6]')
carboxylatePatt = Chem.MolFromSmarts('C(=[OD1])-[OX1]')
carboxylPatt = Chem.MolFromSmarts('C(=[OD1])-[OX2]')

specialCases = ((carboxylatePatt, 97), (esterPatt, 96), (carboxylPatt, 98), (amidePatt, 95),
                (ketonePatt, 94), (aldehydePatt, 93))


def SATISTypes(mol, neighborsToInclude=4):
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
  nAtoms = mol.GetNumAtoms()
  atoms = mol.GetAtoms()
  atomicNums = [atom.GetAtomicNum() for atom in atoms]

  # Collect all atom indices that match one of the special patterns
  specialCaseMatches = []
  for patt, specialCaseIdx in specialCases:
    matches = mol.GetSubstructMatches(patt)
    if matches:
      matches = set(itertools.chain(*matches))
      specialCaseMatches.append((specialCaseIdx, matches))

  codes = [None] * nAtoms
  for i, atom in enumerate(atoms):
    code = [99] * (neighborsToInclude + 1)

    # Atom
    code[0] = min(atom.GetAtomicNum(), 99)

    # Get atomic numbers of connected neighbours and use for code
    otherIndices = [x.GetIdx() for x in atom.GetNeighbors()]
    otherNums = sorted([atomicNums[x] for x in otherIndices] + [1] * atom.GetTotalNumHs())
    if len(otherNums) > neighborsToInclude:
      # Get the last neighborsToInclude elements from otherNums
      otherNums = otherNums[-neighborsToInclude:]
    for j, otherNum in enumerate(otherNums, 1):
      code[j] = min(otherNum, 99)

    # Handle special cases where we have less than neighborsToInclude neighbors
    if len(otherNums) < neighborsToInclude and code[0] in [6, 8]:
      atomIdx = atom.GetIdx()
      for specialCaseIdx, matches in specialCaseMatches:
        if atomIdx in matches:
          code[-1] = specialCaseIdx
          break

    codes[i] = ''.join('%02d' % (x) for x in code)
  return codes
