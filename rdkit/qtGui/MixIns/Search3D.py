#
#  Copyright (C) 2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Enables the Search3D (pharmacophore search) widget

"""
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.Search3D import Search3D

REQUIRED_MIXINS = ['MolBrowser']
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']


METHODS_DEFINED = {
  '__LaunchSearchWidget':'rdkit.qtGui.GuiBase.GuiBase.s3dLaunchSearchWidget',
  }


def __LaunchSearchWidget(self):
  """ launches or displays our _Search3DWidget_
    (method of _GuiBase_)    

  """
  if self._s3dSearchWidget is None:
    self._s3dSearchWidget = Search3D.Search3DApp()
    self._s3dSearchWidget.show() 
    self._s3dSearchWidget.searchWidget.debugInit()
  else:
    self._s3dSearchWidget.show() 
    

def LocalInit(self):
  self._s3dSearchWidget=None
  self._mbMolMenu.insertItem(self.trUtf8('&Pharmacophore Search'),
                             self.s3dLaunchSearchWidget)




                                
                                
