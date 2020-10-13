from rdkit import RDConfig
import unittest
import random
from os import environ
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
import numpy as np


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)

  def test2(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DSVG(300, 300, -1, -1, True)
    do = d.drawOptions()
    do.atomLabels[3] = 'foolabel'
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find(">f</text>") != -1)
    self.assertTrue(txt.find(">o</text>") != -1)
    self.assertTrue(txt.find(">l</text>") != -1)
    self.assertTrue(txt.find(">a</text>") != -1)

  @unittest.skipUnless(hasattr(Draw, 'MolDraw2DCairo'), 'Cairo support not enabled')
  def testGithubIssue571(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DCairo(300, 300)
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()

  def testPrepareForDrawing(self):
    m = Chem.MolFromSmiles('c1ccccc1[C@H](F)Cl')
    nm = rdMolDraw2D.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetNumAtoms(), 9)
    self.assertEqual(nm.GetNumConformers(), 1)
    m = Chem.MolFromSmiles('C1CC[C@H]2NCCCC2C1')
    nm = rdMolDraw2D.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetNumAtoms(), 11)
    self.assertEqual(nm.GetNumConformers(), 1)
    nm = rdMolDraw2D.PrepareMolForDrawing(m, addChiralHs=False)
    self.assertEqual(nm.GetNumAtoms(), 10)
    self.assertEqual(nm.GetNumConformers(), 1)

  def testRepeatedPrepareForDrawingCalls(self):
    m = Chem.MolFromMolBlock("""
          11280715312D 1   1.00000     0.00000     0

 33 36  0     1  0            999 V2000
    7.6125   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
    7.0917   -6.0917    0.0000 C   0  0  1  0  0  0           0  0  0
    6.4792   -6.8917    0.0000 C   0  0  2  0  0  0           0  0  0
    8.1292   -6.0792    0.0000 N   0  0  0  0  0  0           0  0  0
    5.5042   -6.8917    0.0000 C   0  0  3  0  0  0           0  0  0
   11.2375   -4.8542    0.0000 N   0  0  0  0  0  0           0  0  0
    9.6792   -5.1667    0.0000 N   0  0  3  0  0  0           0  0  0
    5.9917   -6.5417    0.0000 C   0  0  0  0  0  0           0  0  0
    7.6042   -5.1917    0.0000 O   0  0  0  0  0  0           0  0  0
   10.7167   -5.1625    0.0000 C   0  0  0  0  0  0           0  0  0
    6.2917   -7.4667    0.0000 C   0  0  0  0  0  0           0  0  0
    6.5750   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
   10.2000   -4.8667    0.0000 C   0  0  0  0  0  0           0  0  0
    8.6500   -5.7792    0.0000 C   0  0  3  0  0  0           0  0  0
    8.6417   -5.1792    0.0000 C   0  0  0  0  0  0           0  0  0
    9.1667   -6.0750    0.0000 C   0  0  0  0  0  0           0  0  0
    9.6875   -5.7667    0.0000 C   0  0  0  0  0  0           0  0  0
    9.1542   -4.8750    0.0000 C   0  0  0  0  0  0           0  0  0
    5.6917   -7.4667    0.0000 C   0  0  0  0  0  0           0  0  0
    5.2042   -7.4042    0.0000 F   0  0  0  0  0  0           0  0  0
    4.9875   -6.5917    0.0000 F   0  0  0  0  0  0           0  0  0
    7.5167   -6.5167    0.0000 O   0  0  0  0  0  0           0  0  0
   11.7542   -5.1500    0.0000 C   0  0  0  0  0  0           0  0  0
   11.2417   -6.0542    0.0000 C   0  0  0  0  0  0           0  0  0
   10.7250   -5.7625    0.0000 C   0  0  0  0  0  0           0  0  0
    6.5750   -5.1917    0.0000 C   0  0  0  0  0  0           0  0  0
    6.0542   -6.0917    0.0000 C   0  0  0  0  0  0           0  0  0
   11.7667   -5.7542    0.0000 C   0  0  0  0  0  0           0  0  0
   12.2750   -4.8417    0.0000 C   0  0  0  0  0  0           0  0  0
    6.0542   -4.8917    0.0000 C   0  0  0  0  0  0           0  0  0
    5.5375   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
    5.5375   -5.1917    0.0000 C   0  0  0  0  0  0           0  0  0
    6.3167   -6.3042    0.0000 H   0  0  0  0  0  0           0  0  0
  2  1  1  0     0  0
  3  2  1  0     0  0
  4  1  1  0     0  0
  5  8  1  0     0  0
  6 10  1  0     0  0
  7 17  1  0     0  0
  8  3  1  0     0  0
  9  1  2  0     0  0
 10 13  1  0     0  0
 11  3  1  0     0  0
 12  2  1  0     0  0
 13  7  1  0     0  0
 14  4  1  0     0  0
 15 14  1  0     0  0
 16 14  1  0     0  0
 17 16  1  0     0  0
 18 15  1  0     0  0
 19 11  1  0     0  0
 20  5  1  0     0  0
 21  5  1  0     0  0
  2 22  1  6     0  0
 23  6  2  0     0  0
 24 25  1  0     0  0
 25 10  2  0     0  0
 26 12  1  0     0  0
 27 12  2  0     0  0
 28 24  2  0     0  0
 29 23  1  0     0  0
 30 26  2  0     0  0
 31 27  1  0     0  0
 32 31  2  0     0  0
  3 33  1  6     0  0
  7 18  1  0     0  0
 19  5  1  0     0  0
 32 30  1  0     0  0
 28 23  1  0     0  0
M  END""")
    nm = Draw.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondDir(), Chem.BondDir.NONE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondDir(), Chem.BondDir.BEGINWEDGE)
    nm = Draw.PrepareMolForDrawing(nm)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondDir(), Chem.BondDir.NONE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondDir(), Chem.BondDir.BEGINWEDGE)

  def testDrawMoleculesArgs(self):
    smis = [
      'O=C1c2cccc3c(N4CCCCC4)ccc(c23)C(=O)N1c1ccccn1', 'Cc1ccc[n+](C2=Nc3ccccc3[N-]C2=C(C#N)C#N)c1',
      'CC(=O)NC1=NN(C(C)=O)C(c2ccccn2)S1', 'COc1cc(Cc2nccc3cc(OC)c(OC)cc23)c(S(=O)(=O)O)cc1OC'
    ]
    tms = [Chem.MolFromSmiles(x) for x in smis]
    [rdMolDraw2D.PrepareMolForDrawing(x) for x in tms]
    drawer = rdMolDraw2D.MolDraw2DSVG(600, 600, 300, 300)
    p = Chem.MolFromSmarts('c1ccccn1')
    matches = [x.GetSubstructMatch(p) for x in tms]
    acolors = []
    for mv in matches:
      clrs = {}
      random.seed(0xf00d)
      for idx in mv:
        clrs[idx] = (random.random(), random.random(), random.random())
      acolors.append(clrs)
    # the bonds between the matching atoms too
    bnds = []
    for mol, mv in zip(tms, matches):
      tmp = []
      for bnd in p.GetBonds():
        tmp.append(
          mol.GetBondBetweenAtoms(mv[bnd.GetBeginAtomIdx()], mv[bnd.GetEndAtomIdx()]).GetIdx())
      bnds.append(tmp)
    drawer.DrawMolecules(tms, highlightAtoms=matches, highlightBonds=bnds,
                         highlightAtomColors=acolors)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # 4 molecules, 6 bonds each:
    self.assertEqual(svg.count('fill:none;fill-rule:evenodd;stroke:#FF7F7F'), 24)
    # 4 molecules, one atom each:
    self.assertEqual(svg.count('fill:#DB2D2B;fill-rule:evenodd;stroke:#DB2D2B'), 4)

  def testGetDrawCoords(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(m)
    conf = m.GetConformer()
    for idx in range(m.GetNumAtoms()):
      pos = conf.GetAtomPosition(idx)
      pos = Geometry.Point2D(pos.x, pos.y)
      dpos1 = d.GetDrawCoords(idx)
      dpos2 = d.GetDrawCoords(pos)
      self.assertAlmostEqual(dpos1.x, dpos2.x, 6)
      self.assertAlmostEqual(dpos1.y, dpos2.y, 6)

  def testReaction1(self):
    rxn = AllChem.ReactionFromSmarts(
      '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',
      useSmiles=True)
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)
    #print(txt,file=open('blah1.svg','w+'))

  def testReaction2(self):
    rxn = AllChem.ReactionFromSmarts(
      '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',
      useSmiles=True)
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn, highlightByReactant=True)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)
    #print(txt,file=open('blah2.svg','w+'))

  def testReaction3(self):
    rxn = AllChem.ReactionFromSmarts(
      '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',
      useSmiles=True)
    colors = [(0.3, 0.7, 0.9), (0.9, 0.7, 0.9), (0.6, 0.9, 0.3), (0.9, 0.9, 0.1)]
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn, highlightByReactant=True, highlightColorsReactants=colors)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)

  def testReaction4(self):
    rxn = AllChem.ReactionFromSmarts(
      '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',
      useSmiles=True)
    colors = [(100, 155, 245), (0, 45, 155)]
    d = Draw.MolDraw2DSVG(900, 300)
    self.assertRaises(ValueError, d.DrawReaction, rxn, True, colors)

  def testBWDrawing(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") >= 0)
    self.assertTrue(txt.find("stroke:#00CC00") >= 0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().useBWAtomPalette()
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") >= 0)
    self.assertTrue(txt.find("stroke:#00CC00") == -1)

  def testUpdatePalette(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") >= 0)
    self.assertTrue(txt.find("stroke:#00CC00") >= 0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().updateAtomPalette({6: (1, 1, 0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") == -1)
    self.assertTrue(txt.find("stroke:#00CC00") >= 0)
    self.assertTrue(txt.find("stroke:#FFFF00") >= 0)

  def testSetPalette(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") >= 0)
    self.assertTrue(txt.find("stroke:#00CC00") >= 0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().setAtomPalette({-1: (1, 1, 0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") == -1)
    self.assertTrue(txt.find("stroke:#00CC00") == -1)
    self.assertTrue(txt.find("stroke:#FFFF00") >= 0)

    # try a palette that doesn't have a default:
    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().setAtomPalette({0: (1, 1, 0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000") >= 0)
    self.assertTrue(txt.find("stroke:#00CC00") == -1)
    self.assertTrue(txt.find("stroke:#FFFF00") == -1)

  def testGithub1829(self):
    d = Draw.MolDraw2DSVG(300, 300, 100, 100)
    d.DrawMolecules(tuple())
    d.FinishDrawing()
    d.GetDrawingText()

  def testSetLineWidth(self):
    " this was github #2149 "
    m = Chem.MolFromSmiles('CC')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke-width:2.0px") >= 0)
    self.assertTrue(txt.find("stroke-width:4.0px") == -1)
    d = Draw.MolDraw2DSVG(300, 300)
    d.SetLineWidth(4)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke-width:2.0px") == -1)
    self.assertTrue(txt.find("stroke-width:4.0px") >= 0)

  def testPrepareAndDrawMolecule(self):
    m = Chem.MolFromSmiles("C1N[C@@H]2OCC12")
    d = Draw.MolDraw2DSVG(300, 300, -1, -1, True)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find(">H</text>") > 0)

  def testAtomTagging(self):
    m = Chem.MolFromSmiles("C1N[C@@H]2OCC12")
    d = Draw.MolDraw2DSVG(300, 300)
    dm = Draw.PrepareMolForDrawing(m)
    rdMolDraw2D.PrepareAndDrawMolecule(d, dm)
    d.TagAtoms(dm, events={'onclick': 'alert'})
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<circle") > 0)
    self.assertTrue(txt.find("onclick=") > 0)

  def testMolContours(self):
    m = Chem.MolFromSmiles("C1N[C@@H]2OCC12")
    dm = Draw.PrepareMolForDrawing(m)

    conf = dm.GetConformer()
    gs = []
    ws = []
    hs = []
    for i in range(conf.GetNumAtoms()):
      p = conf.GetAtomPosition(i)
      p2 = Geometry.Point2D(p.x, p.y)
      gs.append(p2)
      hs.append(0.4)
      if not i % 2:
        ws.append(-0.5)
      else:
        ws.append(1)

    d = Draw.MolDraw2DSVG(300, 300)
    d.ClearDrawing()
    Draw.ContourAndDrawGaussians(d, gs, hs, ws, mol=dm)
    d.drawOptions().clearBackground = False
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    with open("contour_from_py_1.svg", 'w+') as outf:
      print(txt, file=outf)

    d = Draw.MolDraw2DSVG(300, 300)
    d.ClearDrawing()
    ps = Draw.ContourParams()
    ps.fillGrid = True
    Draw.ContourAndDrawGaussians(d, gs, hs, ws, params=ps, mol=dm)
    d.drawOptions().clearBackground = False
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    with open("contour_from_py_2.svg", 'w+') as outf:
      print(txt, file=outf)

  def testGridContours(self):
    grid = np.zeros((50, 100), np.double)
    ycoords = list(np.arange(0, 5, 0.1))
    xcoords = list(np.arange(0, 10, 0.1))
    for i in range(grid.shape[1]):
      gxp = xcoords[i]
      for j in range(grid.shape[0]):
        gyp = ycoords[j]
        for v in range(1, 4):
          dvx = 2 * v - gxp
          dvy = v - gyp
          d2 = dvx * dvx + dvy * dvy
          if d2 > 0:
            grid[j, i] += 1 / np.sqrt(d2)

    # sg = 255 * grid / np.max(grid)
    # from PIL import Image
    # img = Image.fromarray(sg.astype(np.uint8))
    # img.save('img.png')

    d = Draw.MolDraw2DSVG(300, 300)
    d.ClearDrawing()
    Draw.ContourAndDrawGrid(d, np.transpose(grid), xcoords, ycoords)
    d.drawOptions().clearBackground = False
    d.FinishDrawing()
    txt = d.GetDrawingText()
    with open("contour_from_py_3.svg", 'w+') as outf:
      print(txt, file=outf)

  def testExtraDrawingCommands(self):
    " this is another test just to make sure that things work "
    m = Chem.MolFromMolBlock("""
  Mrv1810 07271915232D          

  6  6  0  0  0  0            999 V2000
   -1.5      -1.5       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5       0.0       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0       0.0       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0      -1.5       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5       1.5       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5      -1.5       0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  1  1  0  0  0  0
  2  5  1  0  0  0  0
  4  6  1  0  0  0  0
M  END
""")
    d = Draw.MolDraw2DSVG(300, 300)
    dm = Draw.PrepareMolForDrawing(m)
    d.DrawMolecule(dm)
    conf = dm.GetConformer()
    ps3 = (
      conf.GetAtomPosition(0),
      conf.GetAtomPosition(1),
      conf.GetAtomPosition(3),
    )
    ps = [Geometry.Point2D(p.x, p.y) for p in ps3]
    d.DrawPolygon(ps)
    d.DrawArrow(Geometry.Point2D(0, 0), Geometry.Point2D(0, 1.5))
    d.DrawAttachmentLine(Geometry.Point2D(0, 0), Geometry.Point2D(1.5, 0), (0.5, 0.5, 0.5), len=0.5)
    d.DrawLine(Geometry.Point2D(0, 0), Geometry.Point2D(1.5, 0))
    d.DrawWavyLine(Geometry.Point2D(0, 0), Geometry.Point2D(1.5, 1.5), (0, 0, 0), (1, 0.2, 0.2))
    d.FinishDrawing()
    txt = d.GetDrawingText()
    with open("extras_1.svg", "w+") as outf:
      outf.write(txt)

  def testSetDrawOptions(self):
    m = Chem.MolFromSmiles('CCNC(=O)O')
    d = rdMolDraw2D.MolDraw2DSVG(250, 200, -1, -1, True)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertNotEqual(txt.find("fill:#0000FF' >N</text>"), -1)
    self.assertEqual(txt.find("fill:#000000' >N</text>"), -1)

    d = rdMolDraw2D.MolDraw2DSVG(250, 200, -1, -1, True)
    do = rdMolDraw2D.MolDrawOptions()
    do.useBWAtomPalette()
    d.SetDrawOptions(do)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertEqual(txt.find("fill:#0000FF' >N</text>"), -1)
    self.assertNotEqual(txt.find("fill:#000000' >N</text>"), -1)

  def testAlternativeFreetypeFont(self):
    # this one, you have to look at the pictures
    m = Chem.MolFromSmiles('S(=O)(=O)(O)c1c(Cl)c(Br)c(I)c(F)c(N)1')
    d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    with open('test_ff.svg', 'w') as f:
      f.write(txt)

    d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    do = rdMolDraw2D.MolDrawOptions()
    rdbase = environ['RDBASE']
    if rdbase:
      do.fontFile = '{}/Code/GraphMol/MolDraw2D/Amadeus.ttf'.format(rdbase)
      d.SetDrawOptions(do)
      rdMolDraw2D.PrepareAndDrawMolecule(d, m)
      d.FinishDrawing()
      txt = d.GetDrawingText()
      with open('test_aff.svg', 'w') as f:
        f.write(txt)
    else:
      pass

  def testExplictMethyl(self):
    m = Chem.MolFromSmiles('CC')
    d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertEqual(txt.find("class='atom-"), -1)

    d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    do = rdMolDraw2D.MolDrawOptions()
    do.explicitMethyl = True
    d.SetDrawOptions(do)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertNotEqual(txt.find("class='atom-"), -1)

  def testDrawMoleculeWithHighlights(self):
    COLS = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.55, 0.0)]

    def get_hit_atoms_and_bonds(mol, smt):
      alist = []
      blist = []
      q = Chem.MolFromSmarts(smt)
      for match in mol.GetSubstructMatches(q):
        alist.extend(match)

      for ha1 in alist:
        for ha2 in alist:
          if ha1 > ha2:
            b = mol.GetBondBetweenAtoms(ha1, ha2)
            if b:
              blist.append(b.GetIdx())

      return alist, blist

    def add_colours_to_map(els, cols, col_num):
      for el in els:
        if not el in cols:
          cols[el] = []
        if COLS[col_num] not in cols[el]:
          cols[el].append(COLS[col_num])

    def do_a_picture(smi, smarts, label):

      rdDepictor.SetPreferCoordGen(False)
      mol = Chem.MolFromSmiles(smi)
      mol = Draw.PrepareMolForDrawing(mol)

      acols = {}
      bcols = {}
      h_rads = {}
      h_lw_mult = {}

      for i, smt in enumerate(smarts):
        alist, blist = get_hit_atoms_and_bonds(mol, smt)
        h_rads[alist[0]] = 0.4
        h_lw_mult[blist[0]] = 2
        col = i % 4
        add_colours_to_map(alist, acols, col)
        add_colours_to_map(blist, bcols, col)

      d = rdMolDraw2D.MolDraw2DSVG(500, 500)
      d.drawOptions().fillHighlights = False
      d.DrawMoleculeWithHighlights(mol, label, acols, bcols, h_rads, h_lw_mult, -1)

      d.FinishDrawing()
      return d.GetDrawingText()

    smi = 'CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]'
    smarts = ['CONN', 'N#CC~CO', 'C=CON', 'CONNCN']
    txt = do_a_picture(smi, smarts, 'pyTest2')
    self.assertGreater(txt.find('stroke:#FF8C00;stroke-width:5.5'), -1)
    self.assertEqual(
      txt.find("ellipse cx='244.253' cy='386.518'"
               " rx='11.9872' ry='12.8346'"
               " style='fill:none;stroke:#00FF00'"), -1)

    # test for no-longer-mysterious OSX crash.
    smi = 'c1ccccc1Cl'
    smarts = []
    do_a_picture(smi, smarts, 'pyTest3')

  @unittest.skipUnless(hasattr(Draw, 'MolDraw2DCairo'), 'Cairo support not enabled')
  def testPNGMetadata(self):
    m = Chem.MolFromMolBlock('''
  Mrv2014 08172015242D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.31 -1.3337 0 0
M  V30 2 C 3.6437 -2.1037 0 0
M  V30 3 O 4.9774 -1.3337 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END''')
    d = Draw.MolDraw2DCairo(200, 200)
    d.DrawMolecule(m)
    txt = d.GetDrawingText()
    nm = Chem.MolFromPNGString(txt)
    self.assertEqual(Chem.MolToSmiles(m), Chem.MolToSmiles(nm))


if __name__ == "__main__":
  unittest.main()
