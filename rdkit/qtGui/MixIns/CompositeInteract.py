#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Mixin for interacting with composite models

"""
from rdkit import RDConfig
from qt import *
from rdkit.qtGui import GuiTextViewer
from rdkit.qtGui.GuiLib.CompositeWindowImpl import CompositeWindow
from rdkit.qtGui.GuiLib.DecTreeWindow import TreeWindow
from rdkit.qtGui import qtUtils
import cPickle
import os.path,types

REQUIRED_MIXINS = ['PiddleCanvas','CompositeBase']
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']
METHODS_DEFINED = {
  '__LoadPickledComposite':'rdkit.qtGui.GuiBase.GuiBase.ciLoadPickledComposite',
  '__LoadCompositesFromDb':'rdkit.qtGui.GuiBase.GuiBase.ciLoadCompositesFromDb',
  '__LaunchCompositeBrowser':'rdkit.qtGui.GuiBase.GuiBase.ciLaunchCompositeBrowser',
  '__AddComposite':'rdkit.qtGui.GuiBase.GuiBase.ciAddComposite',
  }
VARS_TO_SAVE = [
  'self._ciComposites',
]

def __AddComposite(self,composite,name=None):
  """ adds a composite to a tab in our compositeWindow
    (method of _GuiBase_)

  **Arguments**

    - composite: the composite to add

    - name: (optional) the name to use for the new tab

  **Notes**

    - this causes our compositeWindow to be displayed (if it is not
      already) 

  """
  if composite:
    self._ciComposites.append(composite)
    self._ciWindow.addComposite(composite,drawTarget=self._ciTreeWin,name=name)
    self._ciWindow.show()
  
                   
def __LoadPickledComposite(self,idx=0,fileN=None):
  """ loads a pickled composite and displays it
    (method of _GuiBase_)    

  **Arguments**

    - idx: (optional) not used

    - fileN: (optional) the filename from which to load the
      composite.  If this is not provided, a _QFileDialog_ will be
      launched to prompt for the name

  **Notes**

    - heavy lifting is done by _self.ciAddComposite()_

      
  """
  if not fileN:
    fileN = str(QFileDialog.getOpenFileName('.','Pickled files (*.pkl);;All files (*.*)'))
  if not fileN:
    return None
  else:
    try:
      inF = open(fileN,'rb')
    except IOError:
      qtUtils.error('could not open file %s for reading'%(fileN))
      return None
    compos = cPickle.load(inF)
    nameBase = os.path.split(fileN)[-1]
    self.ciAddComposite(compos,name=nameBase)
    return compos


def __LoadCompositesFromDb(self,idx=0,fileN=None):
  """ #DOC
    (method of _GuiBase_)    

  **Arguments**

    - idx: (optional) not used

  **Notes**

    - heavy lifting is done by _self.ciAddComposite()_

      
  """
  # FIX: get the stupid directory stuff done properly here
  from rdkit.qtGui.DbQueryDialog import DbQueryDialog

  dlg = DbQueryDialog()
  widget = dlg.dbWidget()
  widget.setWhat('model')
  # FIX: implementation specific
  widget.sqlWhatBox.setEnabled(0)

  res = []
  if dlg.exec_loop() == QDialog.Accepted:
    conn = widget.getConn()
    d = conn.GetData(fields=str(widget.sqlWhat()),
                     where=str(widget.sqlWhere()),
                     join=str(widget.sqlJoin()))
    for i in range(len(d)):
      compos = cPickle.loads(str(d[i][0]))
      self.ciAddComposite(compos,name='Composite %d'%(i+1))
      res.append(compos)
  else:
    print 'cancelled'


  return res

def __LaunchCompositeBrowser(self,idx=0):
  """ #DOC
    (method of _GuiBase_)    

  **Arguments**

    - idx: (optional) not used

  """
  # FIX: get the stupid directory stuff done properly here
  from rdkit.qtGui.GuiLib.ModelBrowserDialogImpl import ModelBrowserDialog

  self._ciModelBrowserDlg = ModelBrowserDialog(self)
  self._ciModelBrowserDlg.show()

def _SetupCompositeWindow(self):
  """ INTERNAL USE ONLY

    creates, initializes, then hides a _CompositeWindow_

  """
  CompositeWindow.fileOpen = self.ciLoadPickledComposite
  self._ciWindow = CompositeWindow()
  self._ciWindow.hide()
  win = self._ciWindow 

  
def LocalInit(self):
  """ Initialization

    - Adds Composite->Load menu item

    - Adds Compoiste->Load->From File menu item

    - Adds View->Compoiste menu item

  """
  self._ciComposites = []
  self._ciLoadItem = QPopupMenu(self)
  self._compositeMenu.insertItem(self.trUtf8('&Load'),self._ciLoadItem)
  self._ciLoadItem.insertItem(self.trUtf8('From &File'),self.ciLoadPickledComposite)
  self._ciLoadItem.insertItem(self.trUtf8('From &Db'),self.ciLoadCompositesFromDb)
  self._ciLoadItem.insertItem(self.trUtf8('Using &Browser'),self.ciLaunchCompositeBrowser)

  _SetupCompositeWindow(self)
  self._viewMenu.insertItem(self.trUtf8('&Composite'),self._ciWindow.show)
  self._ciTreeWin=TreeWindow()

  self._ciExampleListCtrls = []  

  self._ciScreenWin = GuiTextViewer.GuiTextViewer()
  # FIX: too implementation specific
  self._ciScreenWin.textBrowser.mimeSourceFactory().addFilePath(os.getcwd())

  self._ciModelBrowserDialog=None
