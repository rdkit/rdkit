# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#   All Rights Reserved
#
import sys

from qt import *

if sys.platform == 'win32':
  import win32clipboard
  import win32com.client.gencache

  from rdkit.qtGui.qtActiveX import MakeActiveXClass
  try:
    cdxModule = win32com.client.gencache.EnsureModule("{AF2D2DBA-75E4-4123-BC0B-A57BD5C5C5D2}", 0,
                                                      7, 0)
  except Exception:
    raise ImportError("Chemdraw 6.0 or greater does not appear to be installed.")
else:
  raise ImportError("Chemdraw support only available under Windows")


#----------------------------------------------------------------------
class ChemdrawPanel(QWidget):

  def __init__(self, parent=None, name="test", readOnly=0, size=(300, 300)):
    QWidget.__init__(self, parent, name)
    self.cdx = None
    #self.resize(QSize(300,300))
    self.resize(size[0], size[1])
    # Make a new class that derives from the class in the
    # COM module imported above.  This class also derives from QWidget and
    # implements the machinery needed to integrate the two worlds.
    theClass = MakeActiveXClass(cdxModule.ChemDrawCtl, eventObj=self)

    # Create an instance of that class
    self.cdx = theClass(self)
    if readOnly:
      self.cdx.ViewOnly = 1

    # FIX:
    #  This hackery is due to an apparent problem with PyQt: there is
    #  always a gray box about 30 pixels high in the widget we're deriving
    #  from.
    self.offset = 30
    self.label = QLabel(self, "ChemDraw")
    self.label.setText(name)
    self.label.setAlignment(Qt.AlignHCenter)
    fnt = QApplication.font()
    fnt.setPointSize(14)
    self.label.setFont(fnt)

  def pullData(self, fmt='chemical/daylight-smiles'):
    data = self.cdx.GetData(fmt)
    return str(data)

  def setData(self, data, fmt='chemical/daylight-smiles'):
    self.cdx.Objects.Clear()
    res = self.cdx.SetData(fmt, data)
    return res

  def resizeEvent(self, evt):
    sz = evt.size()
    self.label.setGeometry(0, 0, sz.width(), self.offset)
    self.cdx.MoveWindow((0, self.offset, sz.width(), sz.height()), 1)

  def __del__(self):
    if self.cdx:
      self.cdx = None


# demo code
if __name__ == '__main__':
  import sys

  import container
  a = QApplication(sys.argv)
  widg = QMainWindow()
  panel = ChemdrawPanel(widg)
  panel.show()
  widg.setCentralWidget(panel)
  widg.resize(QSize(300, 300))
  panel.setData("c1ccccc1C(=O)O")
  widg.show()
  a.setMainWidget(widg)
  a.exec_loop()
