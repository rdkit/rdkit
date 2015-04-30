# $Id$
#
#  Copyright (C) 2001-2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the SATIS numbers

"""
from __future__ import print_function
import unittest
from rdkit.Chem import *
from rdkit.Chem import SATIS
from rdkit.six.moves import xrange

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')

  def test1(self):
    """ first set of test cases
    """
    data = [('CC(=O)NC',['0601010106', '0606070895', '0806999995', '0701060699',
                         '0601010107']),
            ('O=CC(=O)C',['0806999993', '0601060893', '0606060894', '0806999994',
                          '0601010106']),
            ('C(=O)OCC(=O)O',['0601080896', '0806999996', '0806069996', '0601010608',
                              '0606080898', '0806999998', '0801069998']),
            ('C(=O)[O-]',['0601080897', '0806999997', '0806999997']),
            ('CP(F)(Cl)(Br)(O)',['0601010115', '1508091735', '0915999999',
                                 '1715999999', '3515999999', '0801159999']),
             ]

    for smi,res in data:
      m = MolFromSmiles(smi)
      satis = SATIS.SATISTypes(m)
      assert satis==res,"Bad SATIS for mol %s: %s should have been %s"%(smi,str(satis),
                                                                        str(res))
    
if __name__ == '__main__':
  unittest.main()

