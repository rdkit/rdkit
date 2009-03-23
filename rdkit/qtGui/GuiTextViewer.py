# $Id$
#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" defines a class for viewing rich text

"""
import sys
from qt import *
from rdkit.qtGui.forms.TextViewer import TextViewer as _Form
from rdkit.qtGui.qtUtils  import warning,error

class GuiTextViewer(_Form):
  """ a class for displaying rich text

  **Form:** _forms.TextViewer_

  """
  def __init__(self,parent=None,name=None):
    _Form.__init__(self,parent,name)
    self._fName=None
    self._dir=None

  def saveIt(self,fName=None,dir=None):
    """ saves the contents of the widget to a file

    **Arguments**

      - fName: (optional) name of the file to use, if this is not
        provided, our own _fName_ attribute will be used, if that is
        not defined, we'll call _self.saveAs_

      - dir: (optional) name of the directory to use, if not provided,
        our own _dir_ attribute will be used.

    """
    import os.path
    if fName is None:
      fName = self._fName
    if dir is None:
      dir = self._dir
      
    if not fName:
      self.saveAs(dir=dir)
    else:
      if not dir:
        dir = '.'
      txt = str(self.textBrowser.text())
      try:
        fullName = os.path.join(dir,fName)
        outF = open(fullName,'w+')
      except IOError:
        error('Problems opening output file %s'%(fullName))
      else:
        outF.write(txt)
        outF.close()

  def saveAs(self,dir=None):
    """ prompts for a filename and then saves the contents of the
        widget to that file 

    **Arguments**

      - dir: (optional) name of the directory to use, if not provided,
        our own _dir_ attribute will be used., if that is not set,
        default to '.'

    **Notes**

      - we get a filename using a _QFileDialog_
      
    """
    import os.path
    if dir is None:
      dir = self._dir
    if not dir:
      dir = '.'

    fileN = str(QFileDialog.getSaveFileName(dir,"HTML files (*.html *.htm);;All files (*.*)"))
    if fileN:
      self._dir,self._fName = os.path.split(fileN)
      self.saveIt(fName=self._fName,dir=self._dir)

      
  def setSource(self,src):
    """ sets the source URL for our internal text browser

    """
    self.textBrowser.setSource(src)
    
  def setText(self,text):
    """ sets the source text for our internal text browser

    """
    self.textBrowser.setText(text)
    

  #
  #  SLOTS
  #
  def fileSave(self):
    """ callback for File->Save

    """
    self.saveIt()
  def fileSaveAs(self):
    """ callback for File->SaveAs

    """
    self.saveAs()

  def editCopy(self):
    """ callback for Edit->Copy

    """
    txt = self.textBrowser.selectedText()
    if txt:
      clip = qApp.clipboard()
      clip.setText(txt)

