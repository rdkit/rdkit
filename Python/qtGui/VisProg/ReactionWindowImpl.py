#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for ReactionWindows (for editing reaction
    node properties)

"""    
from utils import PilTools
from qt import *
from qtGui.VisProg.forms.ReactionWindow import ReactionWindow as _Form
import VPUtils


class ReactionWindow(_Form):
  def __init__(self,rxnObj,rxnFile=None,rxnImg=None):
    """

      Arguments:
        - rxnFile: (optional) name of the reaction file (CDXML at the moment)
        - rxnImg:  (optional) a Pil image

    """
    _Form.__init__(self)

    self.guts = rxnObj
    
    # EFF: it's dumb to keep reparsing this file
    self.rxnFile = rxnFile
    self._nReact,self._nProd=VPUtils.quickCDXMLParse(self.rxnFile)

    self.rxnImg = rxnImg
    if rxnImg is not None:
      self._initPicture(rxnImg)
    self._initInfoFrame()
    self._items = []
    
  def _initPicture(self,rxnImg):
    """ resizes the image, converts it to a pixmap, and inserts it in our
      picture box

    """
    picture = self.rxnPicture
    sz = picture.size().width(),picture.size().height()
    img = PilTools.FitImage(rxnImg,sz)
    self.pm = PilTools.PilImgToQPixmap(img)
    picture.setPixmap(self.pm)

  def _initInfoFrame(self):
    """ sets up the extra information displayed about the reaction

    """
    return
    frame = self.infoFrame
    nCols = self._nReact+self._nProd
    #self._infoLayout = QHBoxLayout(frame)
    self._infoLayout = QGridLayout(frame,1,nCols)

    nSoFar=0
    for i in range(self._nReact):
      lab = QLabel(self.trUtf8('Reactant %d'%(i+1)),frame)
      self._infoLayout.addWidget(lab,0,nSoFar)
      lab = QLabel(self.trUtf8('Useful Information'),frame)
      self._infoLayout.addWidget(lab,1,nSoFar)
      nSoFar+=1

    for i in range(self._nProd):
      lab = QLabel(self.trUtf8('Product %d'%(i+1)),frame)
      self._infoLayout.addWidget(lab,0,nSoFar)
      lab = QLabel(self.trUtf8('Useful Information'),frame)
      self._infoLayout.addWidget(lab,1,nSoFar)
      nSoFar+=1

  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def rxnClick(self):
    """  click on reaction picture """
    print 'clickity click'
    
  def multiStepClick(self):
    """  click on multistep checkbox """
    print 'click'
    enVal = not self.multiProdCheck.isChecked()
    #self.maxCyclesLabel.setEnabled(enVal)
    #self.maxCyclesSpin.setEnabled(enVal)
    #self.intermediatesCheck.setEnabled(enVal)
      
  def accept(self):
    """ accept changes made to the reaction """
    txt = self.nodeTextBox.text()
    self.guts.setLabelText(txt)
    # FIX: uh, bit too specific here
    self.guts._canvas.update()
    _Form.accept(self)
