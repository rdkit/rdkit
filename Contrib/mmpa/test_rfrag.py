#
#  Copyright (C) 2017 greg Landrum
#
#   @@ All Rights Reserved  @@
#
import os
import subprocess
import unittest

from rdkit import RDConfig


class TestCase(unittest.TestCase):

  def test1Github1406(self):
    with open('data/simple.smi') as inf:
      p = subprocess.run(('python', 'rfrag.py'), stdin=inf, stdout=subprocess.PIPE)
    self.assertFalse(p.returncode)
    self.assertEqual(p.stdout, b'''c1ccccc1,benzene,,
Cc1ccccc1,toluene,,C[*:1].c1ccc(cc1)[*:1]
''')


if __name__ == '__main__':
  import sys
  if sys.hexversion >= 0x3050000:
    unittest.main()
  else:
    print('Python >=3.5 required to run tests')
