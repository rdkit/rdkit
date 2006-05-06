#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Allows calculation of molecular descriptors

"""
import RDConfig
from qt import *
from qtGui import GuiTextViewer
from qtGui.GuiLib.MolDescriptorWinImpl import MolDescriptorWin
import os

REQUIRED_MIXINS = ['PiddleCanvas','MolBrowser']
MODULES_ALTERED = ['qtGui.GuiBase']
METHODS_DEFINED = {
  '__LaunchDescriptorWin':'qtGui.GuiBase.GuiBase.mdLaunchDescriptorWin',
  }

def __LaunchDescriptorWin(self):
  if self._mdDescWin is None:
    dlg = MolDescriptorWin(parent=self)
    self._mdDescWin = dlg
  else:
    dlg = self._mdDescWin
  dlg.setDrawTarget(self._pcCanvasWin)
  dlg.show()
  
def LocalInit(self):
  self._mdDescWin = None

  self._mbMolMenu.insertItem(self.trUtf8('&Descriptors'),self.mdLaunchDescriptorWin)


                                
                                
