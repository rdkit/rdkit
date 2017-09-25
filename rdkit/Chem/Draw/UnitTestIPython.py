#
#  Copyright (C) 2016  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for IPython/Jupyter integration
"""
import unittest
from rdkit import Chem
from rdkit.Chem import Draw
try:
  from rdkit.Chem.Draw import IPythonConsole
  from IPython.core.display import SVG
except ImportError:
  IPythonConsole = None


class TestCase(unittest.TestCase):

  def setUp(self):
    if IPythonConsole is not None and Draw.MolsToGridImage != IPythonConsole.ShowMols:
      IPythonConsole.InstallIPythonRenderer()
    self.mol = Chem.MolFromSmiles('c1c(C[15NH3+])ccnc1[C@](Cl)(Br)[C@](Cl)(Br)F')

  def tearDown(self):
    if IPythonConsole is not None:
      IPythonConsole.UninstallIPythonRenderer()

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testGithub1089(self):
    m = Chem.MolFromSmiles('CCCC')

    IPythonConsole.ipython_useSVG = False
    res = Draw.MolsToGridImage((m, m))
    self.assertFalse(isinstance(res, SVG))

    res = Draw.MolsToGridImage((m, m), useSVG=True)
    self.assertTrue(isinstance(res, SVG))

    IPythonConsole.ipython_useSVG = True
    res = Draw.MolsToGridImage((m, m))
    self.assertTrue(isinstance(res, SVG))


if __name__ == '__main__':
  unittest.main()
