# $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Mixin for interacting with cluster trees

#DOC entire file

"""
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.GuiLib.BuildClusterDialogImpl import BuildClusterDialog
from rdkit.qtGui.GuiLib.ClusterWindow import ClusterWindow
from rdkit.qtGui import qtUtils
import cPickle
import os.path,types

REQUIRED_MIXINS = []
MODULES_ALTERED = ['rdkit.qtGui.GuiBase']

METHODS_DEFINED = {
  '__LoadPickledClusterTree':'rdkit.qtGui.GuiBase.GuiBase.clusterLoadPickledTree',
  '__LaunchBuilder':'rdkit.qtGui.GuiBase.GuiBase.clusterLaunchBuilder',
  }
VARS_TO_SAVE = [
]

def __LaunchBuilder(self):
  if not self._clusterBuildDlg:
    self._clusterBuildDlg = BuildClusterDialog()

  self._clusterBuildDlg.show()
  self._clusterBuildDlg.exec_loop()
  
def __LoadPickledClusterTree(self,idx=0,fileN=None):
  w = ClusterWindow()
  w.show()
  clust = w.loadPickledClusterTree(fileN=fileN)
  if clust:
    self._clusterWindows.append(w)
    w.raiseW()
    self._clusterClusters.append(clust)
  else:
    w.hide()
    w.destroy()
    w = None
  
  
def LocalInit(self):
  """ Initialization

    - Adds Cluster->Load menu item

    - Adds Cluster->Load->From File menu item

    - Adds View->Cluster menu item

  """
  self._cluster = None
  self._clusterClusters = []

  self._clusterBuildDlg = None
  
  self._clusterMenu = QPopupMenu(self)
  self.menuBar().insertItem(self.trUtf8('C&lusters'),self._clusterMenu)
  self._clusterLoadItem = QPopupMenu(self)
  self._clusterMenu.insertItem(self.trUtf8('&Load'),self._clusterLoadItem)
  self._clusterLoadItem.insertItem(self.trUtf8('From &File'),self.clusterLoadPickledTree)
  self._clusterMenu.insertItem(self.trUtf8('&Build'),self.clusterLaunchBuilder)

  self._clusterWindows=[]


  self._clusterExampleListCtrls = []  


