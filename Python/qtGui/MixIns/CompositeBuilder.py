#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Mixin for building composite models

"""
import RDConfig
from qt import *
from qtGui.GuiLib.BuildCompositeDialogImpl import BuildCompositeDialog
from ML import BuildComposite

REQUIRED_MIXINS = ['CompositeBase','CompositeInteract']
MODULES_ALTERED = ['qtGui.GuiBase','qtGui.GuiLib.BuildCompositeDialogImpl']
METHODS_DEFINED = {
  '__LaunchBuilder':'qtGui.GuiBase.GuiBase.cbLaunchCompositeBuilder',
  '__BuildCompositeFromDetails':'qtGui.GuiBase.GuiBase.cbBuildCompositeFromDetails',
  '__TransferComposite':'qtGui.GuiLib.BuildCompositeDialogImpl.BuildCompositeDialog.inspectIt',
  }
VARS_TO_SAVE = [
  'self._cbComposites',
  'self._cbDetails',
  ]

def __TransferComposite(self):
  """ transfers our current composite to our parent
    (method of _BuildCompositeDialog_)    

    uses our parent's _ciAddComposite_ method
    
  """
  for compos in self._composites:
    self.parent.ciAddComposite(compos)

def __LaunchBuilder(self):
  """ launches or displays our _BuildCompositeDialog_
    (method of _GuiBase_)    

  """
  if self._cbBuilderDlg is None:
    dlg = BuildCompositeDialog(parent=self)
    self._cbBuilderDlg = dlg
  else:
    dlg = self._cbBuilderDlg    
  dlg.show()
  
def __BuildCompositeFromDetails(self,details,silent=0):
  """ constructs a composite from a details object
    (method of _GuiBase_)    

  **Arguments**

    - details: a _ML.CompositeRun.CompositeRun_ instance

    - silent: (optional) if this is zero (the default), a
      _QProgressDialog_ will be displayed to show the progress of the
      build 

  **Returns**

    the composite model constructed

  **Notes**

    this constructs a QProgressDialog to show the progress of the build

  """
  cb = None
  if not silent:
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    dlg = QProgressDialog('Build Progress','Cancel',details.nModels*details.nRuns)
    dlg.setLabelText("Build Progress")
  self._cbDetails = details
  self._cbComposites = []
  try:
    for i in range(details.nRuns):
      cb = lambda x,y=dlg,z=i*details.nModels:y.setProgress(x+z)
      self._cbComposites.append(BuildComposite.RunIt(details,progressCallback=cb,saveIt=0))
  except:
    import traceback
    traceback.print_exc()
    res = []
  else:
    res = self._cbComposites
  if not silent:
    QApplication.restoreOverrideCursor()

  return res

def LocalInit(self):
  """ initialization

    adds Composite->Build menu item

  """
  self._cbBuilderDlg = None
  self._cbComposites = []
  self._cbDetails = None

  self._compositeMenu.insertItem(self.trUtf8('&Build'),self.cbLaunchCompositeBuilder)


  
