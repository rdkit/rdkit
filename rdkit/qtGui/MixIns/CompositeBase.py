#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Base composite interaction mixin

  Doesn't really do much other than add a _Composite_ menu

"""
from rdkit import RDConfig
from qt import *

REQUIRED_MIXINS = []
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']
METHODS_DEFINED = {
  }

def LocalInit(self):
  """ initialization

    adds Composite menu item

  """
  self._composites = []
  self._compositeMenu = QPopupMenu(self)
  self.menuBar().insertItem(self.trUtf8('&Composite'),self._compositeMenu)


