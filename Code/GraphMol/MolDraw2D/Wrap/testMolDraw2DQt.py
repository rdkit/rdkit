from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

import sys
from PyQt5.Qt import *
import sip


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testBasics(self):
    qimg = QImage(250, 200, QImage.Format_RGB32)
    qptr = QPainter(qimg)
    p = sip.unwrapinstance(qptr)
    print("SIP:", p)
    m = Chem.MolFromSmiles('c1ccccc1O')
    Draw.PrepareMolForDrawing(m)
    print('construct')
    d2d = Draw.MolDraw2DFromQPainter(250, 200, p, -1, -1)
    print('done')

    print('DRAW!')
    d2d.DrawMolecule(m)

    qimg.save("foo.png")
    qptr = None
    print('DINE')


if __name__ == "__main__":
  app = QApplication(sys.argv)
  unittest.main()
