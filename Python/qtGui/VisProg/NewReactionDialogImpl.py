#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for NewReactionDialogs

"""    
from qt import *
from qtGui.VisProg.forms.NewReactionDialog import NewReactionDialog as _Form
from qtGui import qtUtils
import os,os.path
import VPLib


class NewReactionDialog(_Form):
  def __init__(self,initDir='',canvas=None):
    _Form.__init__(self)
    self._rxnDir = initDir
    self._imgDir = ''
    self._canvas = canvas
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def spinClick(self):
    """ change one of the spin buttons

      launch a filebrowser for the reaction file name

    """
    self.radio.setButton(2)

  def rxnFileClick(self):
    """ click on browse button for reaction file

      launch a filebrowser for the reaction file name

    """
    self.radio.setButton(0)
    fileN = str(QFileDialog.getOpenFileName(self._rxnDir,"Reactions (*.cdxml)"))
    if fileN:
      self._rxnDir,fileN = os.path.split(fileN)
      self.rxnField.setText(fileN)
      if str(self.imgField.text())=='':
        base = fileN[0:fileN.rfind('.')]
        imgFileN = '%s.gif'%(base)
        self.imgField.setText(imgFileN)
        self._imgDir = self._rxnDir

  def imgFileClick(self):
    """ click on browse button for image file

       launch a filebrowser for the image file name

    """
    fileN = str(QFileDialog.getOpenFileName(self._rxnDir,"Images (*.gif *.jpg *.jpeg)"))
    if fileN:
      self._imgDir,fileN = os.path.split(fileN)
      self.imgField.setText(fileN)

  def accept(self):
    """ close dialog with the ok button

      constructs a new ReactionItem and stores it as self._node
      
    """
    print 'accept'
    self._node = None
    if self.fileRadio.isOn():
      print 'file'
      rxnName = os.path.join(self._rxnDir,str(self.rxnField.text()))
      imgName = os.path.join(self._imgDir,str(self.imgField.text()))
      if rxnName != '':
        self._node = VPLib.VPReactionItem(self._canvas,
                                          rxnFile=rxnName,
                                          rxnImg=imgName)
    elif self.smirksRadio.isOn():
      print 'smirks'
    elif self.numberRadio.isOn():
      print 'number'
      nReact = self.nReactSpinner.value()
      nProd = self.nProdSpinner.value()
      self._node = VPLib.VPReactionItem(self._canvas,
                                    nReactants=nReact,nProducts=nProd)
      print 'NODE:',self._node
                                    
    _Form.accept(self)

