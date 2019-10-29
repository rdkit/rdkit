#
#  Copyright (C) 2001-2017  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


import unittest

from rdkit import Chem
from rdkit.Chem import SATIS


class TestCase(unittest.TestCase):

  def test1(self):
    data = [('CC(=O)NC', ['0601010106', '0606070895', '0806999995', '0701060699', '0601010107']),
            ('O=CC(=O)C', ['0806999993', '0601060893', '0606060894', '0806999994', '0601010106']),
            ('C(=O)OCC(=O)O', ['0601080896', '0806999996', '0806069996', '0601010608', '0606080898',
                               '0806999998', '0801069998']),
            ('C(=O)[O-]', ['0601080897', '0806999997', '0806999997']),
            ('CP(F)(Cl)(Br)(O)', ['0601010115', '1508091735', '0915999999', '1715999999',
                                  '3515999999', '0801159999']),
            ('C=C', ['0601010699', '0601010699']), ]
    msg = "Bad SATIS for mol %s: %s should have been %s"
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      satis = SATIS.SATISTypes(m)
      self.assertEqual(satis, res, msg=msg % (smi, str(satis), str(res)))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
