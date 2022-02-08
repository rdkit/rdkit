#
#  Copyright (C) 2020  Greg Landrum and T5 Informatics GmbH
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import rdMolEnumerator
import unittest


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testLinkNodes(self):
    mb = '''one linknode
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 LINKNODE 1 4 2 1 2 1 5
M  V30 END CTAB
M  END'''
    m = Chem.MolFromMolBlock(mb)
    ps = rdMolEnumerator.MolEnumeratorParams(rdMolEnumerator.EnumeratorType.LinkNode)
    bndl = rdMolEnumerator.Enumerate(m, ps)
    self.assertEqual(bndl.Size(), 4)

  def testPositionVariation(self):
    mb = '''two position variation bonds
  Mrv2007 06242006032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7083 2.415 0 0
M  V30 2 C -3.042 1.645 0 0
M  V30 3 C -3.042 0.105 0 0
M  V30 4 N -1.7083 -0.665 0 0
M  V30 5 C -0.3747 0.105 0 0
M  V30 6 C -0.3747 1.645 0 0
M  V30 7 * -3.042 0.875 0 0
M  V30 8 F -5.0434 0.875 0 0
M  V30 9 * -1.0415 2.03 0 0
M  V30 10 Cl -1.0415 4.34 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(2 2 3) ATTACH=ANY
M  V30 8 1 9 10 ENDPTS=(2 1 6) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
'''
    m = Chem.MolFromMolBlock(mb)
    ps = rdMolEnumerator.MolEnumeratorParams(rdMolEnumerator.EnumeratorType.PositionVariation)
    bndl = rdMolEnumerator.Enumerate(m, ps)
    self.assertEqual(bndl.Size(), 4)

  def testCombined(self):
    mb = '''
  Mrv2014 12212013392D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6 4.7484 0 0
M  V30 2 C -7.3337 3.9784 0 0
M  V30 3 N -7.3337 2.4383 0 0
M  V30 4 C -6 1.6683 0 0
M  V30 5 C -4.6663 2.4383 0 0
M  V30 6 C -4.6663 3.9784 0 0
M  V30 7 C -5.8773 0.0617 0 0
M  V30 8 C -3.2136 0.2013 0 0
M  V30 9 C -4.5052 -0.6374 0 0
M  V30 10 O -4.4246 -2.1753 0 0
M  V30 11 * -6 4.235 0 0
M  V30 12 C -4.845 6.2355 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 8 9
M  V30 8 1 7 9
M  V30 9 1 5 8
M  V30 10 1 7 4
M  V30 11 1 9 10
M  V30 12 1 11 12 ENDPTS=(3 1 2 6) ATTACH=ANY
M  V30 END BOND
M  V30 LINKNODE 1 3 2 9 7 9 8
M  V30 END CTAB
M  END
'''
    m = Chem.MolFromMolBlock(mb)
    bndl = rdMolEnumerator.Enumerate(m)
    self.assertEqual(bndl.Size(), 9)


if __name__ == '__main__':
  unittest.main()
