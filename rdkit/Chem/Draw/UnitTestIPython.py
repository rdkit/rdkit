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
import warnings
from unittest import mock

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
  def testCairolessRenderingFallsBackToSVG(self):
    old_use_svg = IPythonConsole.ipython_useSVG
    old_show_properties = IPythonConsole.ipython_showProperties
    try:
      IPythonConsole.ipython_useSVG = False
      IPythonConsole.ipython_showProperties = True
      self.mol.SetProp('publicprop', 'ppropval')
      with mock.patch.object(IPythonConsole, '_hasCairoRenderer', return_value=False):
        self.assertTrue(IPythonConsole._useSVG())
        with mock.patch.object(IPythonConsole.Draw, '_moltoimg', return_value=b'png'):
          self.assertEqual(IPythonConsole._toPNG(self.mol), b'png')
        self.assertIsNotNone(IPythonConsole._toSVG(self.mol))
        html = IPythonConsole._toHTML(self.mol)
      self.assertIn('<svg', html)
      self.assertIn('publicprop', html)
    finally:
      IPythonConsole.ipython_useSVG = old_use_svg
      IPythonConsole.ipython_showProperties = old_show_properties

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testExplicitPNGWithoutCairoFallsBackToSVGQuietly(self):
    bitInfo = {}
    Chem.RDKFingerprint(self.mol, bitInfo=bitInfo)
    bitId = next(iter(bitInfo))
    drawingCalls = (
      (Draw.MolsToGridImage, ([self.mol],)),
      (Draw.DrawRDKitBit, (self.mol, bitId, bitInfo)),
      (Draw.DrawRDKitBits, ([(self.mol, bitId, bitInfo)],)),
    )
    with mock.patch.object(IPythonConsole, '_hasCairoRenderer', return_value=False):
      for fn, args in drawingCalls:
        with self.subTest(function=fn.__name__):
          with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.assertIsInstance(fn(*args, useSVG=False), SVG)
            self.assertIsInstance(fn(*args), SVG)
            self.assertIsInstance(fn(*args, useSVG=True), SVG)
          self.assertEqual(caught, [])
      with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        self.assertIsNone(IPythonConsole._toPNGRepr(self.mol))
        self.assertIsNone(IPythonConsole._toReactionPNGRepr(None))
        self.assertIsNone(IPythonConsole._toMolBundlePNGRepr(None))
      self.assertEqual(caught, [])

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

  @unittest.skipIf(IPythonConsole is None, 'IPython not available')
  def testSimpleHighlighting(self):
    m = Chem.MolFromSmiles('c1ccccc1O')
    m.__sssAtoms = [0, 1, 2]
    # These two calls should fail if the methods haven't ben updated
    if hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
      png = IPythonConsole._toPNG(m)
    svg = IPythonConsole._toSVG(m)


@unittest.skipIf(IPythonConsole is None, 'IPython not available')
class Test3DDrawing(unittest.TestCase):

  def _makeMolWithConformers(self, smiles='CC', ids=(0, 1)):
    mol = Chem.MolFromSmiles(smiles)
    for offset, confId in enumerate(ids):
      conf = Chem.Conformer(mol.GetNumAtoms())
      conf.SetId(confId)
      for idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(idx, (1.5 * idx, 10.0 * offset, 0.0))
      mol.AddConformer(conf, assignId=False)
    return mol

  def testDrawMols3DStyles(self):
    cases = [(None, ['stick'] * 7), ('sphere', ['sphere'] * 7),
             (['sphere', 'stick'], ['sphere', 'stick']),
             (('stick', 'sphere'), ['stick', 'sphere'])]
    for drawAs, expected in cases:
      with self.subTest(drawAs=drawAs):
        view = mock.Mock()
        mols = [Chem.MolFromSmiles('CC') for _ in expected]
        colors = ('cyanCarbon', 'redCarbon')
        result = IPythonConsole.drawMols3D(mols, view=view, drawAs=drawAs, colors=colors)
        self.assertIs(result, view.show.return_value)
        self.assertEqual(view.addModel.call_count, len(mols))
        calls = []
        for idx, styleName in enumerate(expected):
          style = {'colorscheme': colors[idx % len(colors)]}
          if styleName == 'sphere':
            style['radius'] = 0.5
          calls.append(mock.call({'model': idx}, {styleName: style}))
        self.assertEqual(view.setStyle.call_args_list[-len(mols):], calls)

  def testAddMolToViewConformer(self):
    for ids in [(0, 1), (3, 7)]:
      mol = self._makeMolWithConformers(ids=ids)
      for style in ['cartoon', 'stick']:
        with self.subTest(ids=ids, style=style):
          view = mock.Mock()
          IPythonConsole.addMolToView(mol, view, confId=ids[1], drawAs=style)
          if style == 'cartoon':
            expected = Chem.MolToPDBBlock(mol, confId=ids[1], flavor=0x20 | 0x10)
            view.addModel.assert_called_once_with(expected, 'pdb')
          else:
            expected = Chem.MolToMolBlock(mol, confId=ids[1])
            view.addModel.assert_called_once_with(expected, 'sdf')

  def testAddLargeMolToViewConformer(self):
    mol = self._makeMolWithConformers(smiles='C' * 999, ids=(3, 7))
    view = mock.Mock()
    IPythonConsole.addMolToView(mol, view, confId=7)
    expected = Chem.MolToPDBBlock(mol, confId=7, flavor=0x20 | 0x10)
    view.addModel.assert_called_once_with(expected, 'pdb')

  def testDrawMols3DCartoonUsesPDB(self):
    mol = self._makeMolWithConformers()
    view = mock.Mock()
    IPythonConsole.drawMols3D([mol, mol], view=view, confIds=[1, 0],
                             drawAs=['cartoon', 'stick'])
    expectedPDB = Chem.MolToPDBBlock(mol, confId=1, flavor=0x20 | 0x10)
    expectedSDF = Chem.MolToMolBlock(mol, confId=0)
    self.assertEqual(view.addModel.call_args_list,
                     [mock.call(expectedPDB, 'pdb'), mock.call(expectedSDF, 'sdf')])

  def testDrawMols3DCartoonDefaultConformer(self):
    for ids in [(0, 1), (7, 3), ()]:
      with self.subTest(ids=ids):
        mol = self._makeMolWithConformers(ids=ids)
        view = mock.Mock()
        IPythonConsole.drawMols3D([mol], view=view, drawAs=['cartoon'])
        expected = Chem.MolToPDBBlock(mol, confId=ids[0] if ids else -1,
                                      flavor=0x20 | 0x10)
        view.addModel.assert_called_once_with(expected, 'pdb')


if __name__ == '__main__':
  unittest.main()
