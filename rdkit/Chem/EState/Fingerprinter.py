# # Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""  EState fingerprinting

"""

import numpy

from rdkit.Chem.EState import AtomTypes, EStateIndices


def FingerprintMol(mol):
  """ generates the EState fingerprints for the molecule

  Concept from the paper: Hall and Kier JCICS _35_ 1039-1045 (1995)

  two numeric arrays are returned:
    The first (of ints) contains the number of times each possible atom type is hit
    The second (of floats) contains the sum of the EState indices for atoms of
      each type.

  """
  if AtomTypes.esPatterns is None:
    AtomTypes.BuildPatts()
  esIndices = EStateIndices(mol)

  nPatts = len(AtomTypes.esPatterns)
  counts = numpy.zeros(nPatts, dtype=numpy.int64)
  sums = numpy.zeros(nPatts, dtype=numpy.float64)
  for i, (_, pattern) in enumerate(AtomTypes.esPatterns):
    matches = mol.GetSubstructMatches(pattern, uniquify=1)
    counts[i] = len(matches)
    for match in matches:
      sums[i] += esIndices[match[0]]
  return counts, sums


def _exampleCode():
  """ Example code for calculating E-state fingerprints """
  from rdkit import Chem
  smis = ['CC', 'CCC', 'c1[nH]cnc1CC(N)C(O)=O', 'NCCc1ccc(O)c(O)c1']
  for smi in smis:
    m = Chem.MolFromSmiles(smi)
    print(smi, Chem.MolToSmiles(m))
    types = AtomTypes.TypeAtoms(m)
    for i in range(m.GetNumAtoms()):
      print('%d %4s: %s' % (i + 1, m.GetAtomWithIdx(i).GetSymbol(), str(types[i])))
    es = EStateIndices(m)
    counts, sums = FingerprintMol(m)
    for i in range(len(AtomTypes.esPatterns)):
      if counts[i]:
        name, _ = AtomTypes.esPatterns[i]
        print('%6s, % 2d, % 5.4f' % (name, counts[i], sums[i]))
    for i in range(len(es)):
      print('% 2d, % 5.4f' % (i + 1, es[i]))
    print('--------')


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
