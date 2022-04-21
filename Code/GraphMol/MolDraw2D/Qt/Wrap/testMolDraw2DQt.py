import sys
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2DQt

try:
  if rdMolDraw2DQt.rdkitQtVersion.startswith('6'):
    from PyQt6.QtGui import *
    Format_RGB32 = QImage.Format.Format_RGB32
  else:
    from PyQt5.Qt import *
    Format_RGB32 = QImage.Format_RGB32
except ImportError:
  # If we can't find Qt, there's nothing we can do
  sip = None
else:
  try:
    # Prefer the PyQt-bundled sip
    if rdMolDraw2DQt.rdkitQtVersion.startswith('6'):
      from PyQt6 import sip
    else:
      from PyQt5 import sip
  except ImportError:
    # No bundled sip, try the standalone package
    try:
      import sip
    except ImportError:
      # No sip at all
      sip = None


@unittest.skipIf(sip is None, "skipping tests because pyqt is not installed")
class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testSIPBasics(self):
    m = Chem.MolFromSmiles('c1ccccc1O')
    Draw.PrepareMolForDrawing(m)
    qimg = QImage(250, 200, Format_RGB32)
    with QPainter(qimg) as qptr:
      p = sip.unwrapinstance(qptr)
      d2d = rdMolDraw2DQt.MolDraw2DFromQPainter_(250, 200, p)
      d2d.DrawMolecule(m)
    qimg.save("testImageFromPyQt-1.png")

  def testBasics(self):
    m = Chem.MolFromSmiles('c1ccccc1O')
    Draw.PrepareMolForDrawing(m)
    qimg = QImage(250, 200, Format_RGB32)
    with QPainter(qimg) as qptr:
      d2d = Draw.MolDraw2DFromQPainter(qptr)
      d2d.DrawMolecule(m)
    qimg.save("testImageFromPyQt-2.png")

  def testMemory1(self):

    def testfunc():
      m = Chem.MolFromSmiles('c1ccccc1O')
      Draw.PrepareMolForDrawing(m)
      qimg = QImage(250, 200, Format_RGB32)
      with QPainter(qimg) as qptr:
        p = sip.unwrapinstance(qptr)
        d2d = rdMolDraw2DQt.MolDraw2DFromQPainter_(250, 200, p, -1, -1)
      raise ValueError("expected")

    with self.assertRaises(ValueError):
      testfunc()

  def testMemory2(self):

    def testfunc():
      m = Chem.MolFromSmiles('c1ccccc1O')
      Draw.PrepareMolForDrawing(m)
      qimg = QImage(250, 200, Format_RGB32)
      with QPainter(qimg) as qptr:
        d2d = Draw.MolDraw2DFromQPainter(qptr)
      raise ValueError("expected")

    with self.assertRaises(ValueError):
      testfunc()


if __name__ == "__main__":
  if sip is not None:
    app = QGuiApplication(sys.argv)
  unittest.main()
