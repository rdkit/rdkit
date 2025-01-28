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
  from IPython.core.display import SVG

  from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
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

    if hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
      IPythonConsole.ipython_useSVG = False
      res = Draw.MolsToGridImage((m, m))
      self.assertFalse(isinstance(res, SVG))

    res = Draw.MolsToGridImage((m, m), useSVG=True)
    self.assertTrue(isinstance(res, SVG))

    IPythonConsole.ipython_useSVG = True
    res = Draw.MolsToGridImage((m, m))
    self.assertTrue(isinstance(res, SVG))

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testGithub3101(self):
    m = Chem.MolFromSmiles('CCCC')

    IPythonConsole.ipython_useSVG = True
    IPythonConsole.drawOptions.addAtomIndices = True
    res = Draw.MolsToGridImage((m, m))
    self.assertIn('class="note"', res.data)

    dopts = rdMolDraw2D.MolDrawOptions()
    res = Draw.MolsToGridImage((m, m), drawOptions=dopts)
    self.assertNotIn('class="note"', res.data)

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testMolPropertyRendering(self):
    m = Chem.MolFromSmiles('CCCC')
    self.assertTrue(hasattr(m, '_repr_html_'))

    m.SetProp('_Name', 'testm')
    m.SetProp('_hiddenprop', 'hpvalue')
    m.SetProp('computedprop', 'cpvalue', computed=True)
    m.SetDoubleProp('foo', 12.3)
    m.SetProp('publicprop', 'ppropval')

    IPythonConsole.ipython_showProperties = True
    html = m._repr_html_()
    self.assertIn('publicprop', html)
    self.assertNotIn('testm', html)
    self.assertNotIn('_hiddenprop', html)
    self.assertNotIn('computedprop', html)

    IPythonConsole.ipython_showProperties = False
    html = m._repr_html_()
    self.assertNotIn('publicprop', html)
    self.assertNotIn('_hiddenprop', html)
    self.assertNotIn('computedprop', html)

    for i in range(10):
      m.SetIntProp(f'prop-{i}', i)
    IPythonConsole.ipython_showProperties = True
    IPythonConsole.ipython_maxProperties = 5
    html = m._repr_html_()
    self.assertIn('publicprop', html)
    self.assertIn('prop-1', html)
    self.assertNotIn('prop-8', html)

    IPythonConsole.ipython_maxProperties = -1
    html = m._repr_html_()
    self.assertIn('publicprop', html)
    self.assertIn('prop-1', html)
    self.assertIn('prop-8', html)

    IPythonConsole.ipython_maxProperties = 10
    html = m._repr_html_()
    self.assertIn('publicprop', html)
    self.assertIn('prop-1', html)
    self.assertNotIn('prop-8', html)

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testMolsMatrixToGridImage(self):
    # For testMolsMatrixToLinear and testMolsMatrixToLinear (which test MolsMatrixToGridImage and its helper _MolsNestedToLinear)
    s = "NC(C)C(=O)"
    mol = Chem.MolFromSmiles(s)
    # Set up matrix with oligomer count for the molecules
    # Should produce this grid:
    # NC(C)C(=O)
    #                                NC(C)C(=O)NC(C)C(=O)
    # NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)                      NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)
    repeats = [[1], [0, 2], [3, 0, 4]]

    # Create molecule if there are 1 or more oligomers;
    # otherwise, use None for molecule because drawing is distorted if use Chem.MolFromSmiles("")
    molsMatrix = [[Chem.MolFromSmiles(s * count) if count else None for count in row]
                  for row in repeats]

    legendsMatrix = [[str(count) + " unit(s)" for count in row] for row in repeats]
    d = Draw.MolsMatrixToGridImage(molsMatrix, useSVG=True, legendsMatrix=legendsMatrix)
    self.assertIsNotNone(d)


if __name__ == '__main__':
  unittest.main()
