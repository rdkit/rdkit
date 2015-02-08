from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import Draw,AllChem

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
        self.failUnless(txt.find("<svg:svg")!=-1)
        self.failUnless(txt.find("</svg:svg>")!=-1)
            
if __name__=="__main__":
    unittest.main()
