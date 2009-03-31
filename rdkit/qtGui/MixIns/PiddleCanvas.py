#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Allows calculation of molecular descriptors

"""
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.PiddleWindowImpl import PiddleWindow


REQUIRED_MIXINS = []
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']
METHODS_DEFINED = {
  '__LaunchCanvasWin':'rdkit.qtGui.GuiBase.GuiBase.pcLaunchCanvasWin',
  }

def __LaunchCanvasWin(self,*args,**kwargs):
  if self._pcCanvasWin is None:
    win = PiddleWindow(*args,**kwargs)
    win.resizeCanvas((600,600))
    self._pcCanvasWin = win
  else:
    win = self._pcCanvasWin
  win.show()
  

def LocalInit(self):
  self._pcCanvasWin = PiddleWindow()
  self._pcCanvasWin.resizeCanvas((600,600))
  self._pcCanvasWin.hide()

  self._viewMenu.insertItem(self.trUtf8('C&anvas'),self.pcLaunchCanvasWin)

