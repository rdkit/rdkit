#
#  Copyright (C) 2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Ties the MolBrowser and ClusterInteract widgets together

"""
import RDConfig
from qt import *
from qtGui.GuiLib.MolBrowserImpl import MolBrowser
from qtGui.GuiLib.MolTable import MolTable
from qtGui.GuiLib.ClusterWindow import ClusterWindow
import os,weakref

REQUIRED_MIXINS = ['MolBrowser','ClusterInteract']
MODULES_ALTERED = ['qtGui.GuiLib.MolBrowserImpl',
                   'qtGui.GuiLib.ClusterWindow']

METHODS_DEFINED = {
  '__ClusterHighlightMol':'qtGui.GuiLib.MolBrowserImpl.MolBrowser.mcClusterHighlightMol',
  '__LocateMol':'qtGui.GuiLib.ClusterWindow.ClusterWindow.mcLocateMol',
  '__ClusterWindowAssociateBrowser':'qtGui.GuiLib.ClusterWindow.ClusterWindow.mcAssociateBrowser',
  }

def __ClusterHighlightMol(self,tbl,row,col):
  if not hasattr(self,'_mcClusterWindow') or not self._mcClusterWindow():
    return
  # FIX: this should eventually be more flexible:
  txt = self.inDataTable.text(row,1)
  self._mcClusterWindow().ptsLine.setText(txt)
  self._mcClusterWindow().highlightPointsSlot()

def __LocateMol(self,clustWin,cluster,targetCol=1):
  if not hasattr(self,'_mcMolBrowser') or not self._mcMolBrowser():
    return
  nm = cluster.GetName()
  tbl = self._mcMolBrowser().inDataTable
  found = False
  for row in range(tbl.numRows()):
    txt = tbl.text(row,targetCol)
    if txt==nm:
      found=True
      break
  if found:
    tbl.setCurrentCell(row,0)
    tbl.ensureCellVisible(row,targetCol)

def _clusterWinPickCenters(self,nToPick,method=None,colHeading='Picked?',targetCol=1):
  picks = self._origPickCenters(nToPick,method=method)
  if hasattr(self,'_mcMolBrowser') and self._mcMolBrowser():
    tbl = self._mcMolBrowser().inDataTable
    hdr = tbl.horizontalHeader()
    lab = str(hdr.label(tbl.numCols()-1))
    if lab!=colHeading:
      nC = tbl.numCols()
      tbl.insertColumns(nC)
      lab = hdr.setLabel(nC,colHeading)
    colId = tbl.numCols()-1
    for row in range(tbl.numRows()):
      txt = str(tbl.text(row,targetCol))
      if txt in picks:
        tbl.setText(row,colId,"1")
      else:
        tbl.setText(row,colId,"0")
    
      
  
  return picks
  
  
def __ClusterWindowAssociateBrowser(self):
  gui = qApp.mainWidget()
  if hasattr(gui,'_mbBrowserWin') and gui._mbBrowserWin:
    brows = gui._mbBrowserWin
    if hasattr(brows,'_mcClusterWindow') and \
       type(brows._mcClusterWindow)==weakref.ReferenceType:
      brows._mcClusterWindow().clearPointClickCallbacks()
      brows._mcClusterWindow()._mcMolBrowser = None
      
    brows._mcClusterWindow = weakref.ref(self)
    brows.inDataTable.clearChangedCallbacks()
    brows.inDataTable.addChangedCallback(brows.mcClusterHighlightMol)
    self._mcMolBrowser = weakref.ref(brows)
    self.clearPointClickCallbacks()
    self.addPointClickCallback(self.mcLocateMol)
    
def _clusterWinInitMenubar(self):
  self._origInitMenubar()
  self._actionMenu.insertItem(self.trUtf8('&Associate with Browser'),self.mcAssociateBrowser)

def LocalInit(self):
  ClusterWindow._origInitMenubar=ClusterWindow._initMenubar
  ClusterWindow._initMenubar=_clusterWinInitMenubar
  ClusterWindow._origPickCenters = ClusterWindow.pickCenters
  ClusterWindow.pickCenters = _clusterWinPickCenters

                                
                                
