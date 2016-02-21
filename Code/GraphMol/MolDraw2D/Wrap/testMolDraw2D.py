from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import Draw,AllChem
from rdkit.Chem.Draw import rdMolDraw2D

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1(self) :
        m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
        AllChem.Compute2DCoords(m)
        d = Draw.MolDraw2DSVG(300,300)
        d.DrawMolecule(m)
        d.FinishDrawing()
        txt = d.GetDrawingText()
        self.assertTrue(txt.find("<svg:svg")!=-1)
        self.assertTrue(txt.find("</svg:svg>")!=-1)

    def test2(self) :
        m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
        AllChem.Compute2DCoords(m)
        d = Draw.MolDraw2DSVG(300,300)
        do = d.drawOptions()
        do.atomLabels[3]='foolabel'
        d.DrawMolecule(m)
        d.FinishDrawing()
        txt = d.GetDrawingText()
        self.assertTrue(txt.find("foolabel")!=-1)

    def testGithubIssue571(self) :
        if not hasattr(Draw,'MolDraw2DCairo'):
            return
        m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
        AllChem.Compute2DCoords(m)
        d = Draw.MolDraw2DCairo(300,300)
        d.DrawMolecule(m)
        d.FinishDrawing()
        txt = d.GetDrawingText()

    def testPrepareForDrawing(self):
        m = Chem.MolFromSmiles('c1ccccc1[C@H](F)Cl')
        nm = rdMolDraw2D.PrepareMolForDrawing(m)
        self.assertEqual(nm.GetNumAtoms(),10)
        self.assertEqual(nm.GetNumConformers(),1)

if __name__=="__main__":
    unittest.main()
