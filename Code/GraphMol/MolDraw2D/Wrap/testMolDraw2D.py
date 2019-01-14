from rdkit import RDConfig
import unittest
import random
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry


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
    d = Draw.MolDraw2DSVG(300, 300)
    do = d.drawOptions()
    do.atomLabels[3] = 'foolabel'
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("foolabel") != -1)

  def testGithubIssue571(self):
    if not hasattr(Draw, 'MolDraw2DCairo'):
      return
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
    smis = ['O=C1c2cccc3c(N4CCCCC4)ccc(c23)C(=O)N1c1ccccn1', 'Cc1ccc[n+](C2=Nc3ccccc3[N-]C2=C(C#N)C#N)c1', 'CC(=O)NC1=NN(C(C)=O)C(c2ccccn2)S1', 'COc1cc(Cc2nccc3cc(OC)c(OC)cc23)c(S(=O)(=O)O)cc1OC']
    tms = [Chem.MolFromSmiles(x) for x in smis]
    [rdMolDraw2D.PrepareMolForDrawing(x) for x in tms]
    drawer = rdMolDraw2D.MolDraw2DSVG(600,600,300,300)
    p = Chem.MolFromSmarts('c1ccccn1')
    matches = [x.GetSubstructMatch(p) for x in tms]
    acolors = []
    for mv in matches:
        clrs = {}
        random.seed(0xf00d)
        for idx in mv:
            clrs[idx] = (random.random(),random.random(),random.random())
        acolors.append(clrs)
    # the bonds between the matching atoms too
    bnds = []
    for mol,mv in zip(tms,matches):
        tmp = []
        for bnd in p.GetBonds():
            tmp.append(mol.GetBondBetweenAtoms(mv[bnd.GetBeginAtomIdx()],mv[bnd.GetEndAtomIdx()]).GetIdx())
        bnds.append(tmp)
    drawer.DrawMolecules(tms,highlightAtoms=matches,highlightBonds=bnds,highlightAtomColors=acolors)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # 4 molecules, 6 bonds each:
    self.assertEqual(svg.count('fill:none;fill-rule:evenodd;stroke:#FF7F7F'),24)
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
        pos = Geometry.Point2D(pos.x,pos.y)
        dpos1 = d.GetDrawCoords(idx)
        dpos2 = d.GetDrawCoords(pos)
        self.assertAlmostEqual(dpos1.x,dpos2.x,6)
        self.assertAlmostEqual(dpos1.y,dpos2.y,6)

  def testReaction1(self):
    rxn = AllChem.ReactionFromSmarts('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',useSmiles=True)
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)
    #print(txt,file=open('blah1.svg','w+'))

  def testReaction2(self):
    rxn = AllChem.ReactionFromSmarts('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',useSmiles=True)
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn,highlightByReactant=True)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)
    #print(txt,file=open('blah2.svg','w+'))

  def testReaction3(self):
    rxn = AllChem.ReactionFromSmarts('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',useSmiles=True)
    colors=[(0.3,0.7,0.9),(0.9,0.7,0.9),(0.6,0.9,0.3),(0.9,0.9,0.1)]
    d = Draw.MolDraw2DSVG(900, 300)
    d.DrawReaction(rxn,highlightByReactant=True,highlightColorsReactants=colors)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg") != -1)
    self.assertTrue(txt.find("</svg>") != -1)

  def testReaction4(self):
    rxn = AllChem.ReactionFromSmarts('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=[O:3])[NH:6][CH3:5].[OH2:4]',useSmiles=True)
    colors=[(100,155,245),(0,45,155)]
    d = Draw.MolDraw2DSVG(900, 300)
    self.assertRaises(ValueError, d.DrawReaction, rxn, True, colors)

  def testBWDrawing(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")>=0)
    self.assertTrue(txt.find("stroke:#00CC00")>=0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().useBWAtomPalette()
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")>=0)
    self.assertTrue(txt.find("stroke:#00CC00")==-1)

  def testUpdatePalette(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")>=0)
    self.assertTrue(txt.find("stroke:#00CC00")>=0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().updateAtomPalette({6:(1,1,0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")==-1)
    self.assertTrue(txt.find("stroke:#00CC00")>=0)
    self.assertTrue(txt.find("stroke:#FFFF00")>=0)


  def testSetPalette(self):
    m = Chem.MolFromSmiles('CCOCNCCl')
    dm = Draw.PrepareMolForDrawing(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")>=0)
    self.assertTrue(txt.find("stroke:#00CC00")>=0)

    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().setAtomPalette({-1:(1,1,0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")==-1)
    self.assertTrue(txt.find("stroke:#00CC00")==-1)
    self.assertTrue(txt.find("stroke:#FFFF00")>=0)

    # try a palette that doesn't have a default:
    d = Draw.MolDraw2DSVG(300, 300)
    d.drawOptions().setAtomPalette({0:(1,1,0)})
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke:#000000")>=0)
    self.assertTrue(txt.find("stroke:#00CC00")==-1)
    self.assertTrue(txt.find("stroke:#FFFF00")==-1)

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
    self.assertTrue(txt.find("stroke-width:2px") >= 0)
    self.assertTrue(txt.find("stroke-width:4px") == -1)
    d = Draw.MolDraw2DSVG(300, 300)
    d.SetLineWidth(4)
    d.DrawMolecule(dm)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("stroke-width:2px") == -1)
    self.assertTrue(txt.find("stroke-width:4px") >= 0)
 
  def testPrepareAndDrawMolecule(self):
    m = Chem.MolFromSmiles("C1N[C@@H]2OCC12")
    d = Draw.MolDraw2DSVG(300, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(d,m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<tspan>H</tspan>")>0)


if __name__ == "__main__":
  unittest.main()
