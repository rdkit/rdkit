#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" mixin providing reaction-cascade builder functionality on a visual programming
canvas

"""
import RDConfig
from qt import *
from qtcanvas import *
import math
import VPItems,VPLib
from NewReactionDialogImpl import NewReactionDialog

REQUIRED_MIXINS = ['VPCanvas']
MODULES_ALTERED = ['qtGui.GuiBase']
METHODS_DEFINED = {
  '__NewSupplyNode':'qtGui.GuiBase.GuiBase.vpcNewSupplyNode',
  '__NewFilterNode':'qtGui.GuiBase.GuiBase.vpcNewFilterNode',
  '__NewReactionNode':'qtGui.GuiBase.GuiBase.vpcNewReactionNode',
  '__NewSinkNode':'qtGui.GuiBase.GuiBase.vpcNewSinkNode',

  }

def __NewSupplyNode(self):
  node = VPLib.VPSupplyItem(self._vpcView._canvas)
  node.moveBy(self.width()/2,self.height()/2)
  node.show()
  self._vpcView.items.append(node)
  self._vpcView._canvas.update()
  self._vpcView._registerAddNodeUndo(node)
  
def __NewReactionNode(self):
  try:
    startDir = self._vpcRxnDir
  except:
    startDir = ''
  dlg = NewReactionDialog(initDir=startDir,
                          canvas=self._vpcView._canvas)
  res = dlg.exec_loop()
  print 'RES:',res
  if res and dlg._node:
    print 'adding'
    node = dlg._node
    node.moveBy(self.width()/2,self.height()/2)
    node.show()
    self._vpcView.items.append(node)
    self._vpcView._canvas.update()
    self._vpcView._registerAddNodeUndo(node)
  else:
    print 'nope'
  
def __NewFilterNode(self):
  node = VPLib.VPFilterItem(self._vpcView._canvas)
  node.moveBy(self.width()/2,self.height()/2)
  node.show()
  self._vpcView.items.append(node)
  self._vpcView._canvas.update()

def __NewSinkNode(self):
  node = VPLib.VPSinkItem(self._vpcView._canvas)
  node.moveBy(self.width()/2,self.height()/2)
  node.show()
  self._vpcView.items.append(node)
  self._vpcView._canvas.update()
  self._vpcView._registerAddNodeUndo(node)
  
def LocalInit(self):
  self._vpcNodeMenu = QPopupMenu(self)
  self.menuBar().insertItem(self.trUtf8('&Nodes'),self._vpcNodeMenu)
  self._vpcNodeMenu.insertItem(self.trUtf8('&Supply'),self.vpcNewSupplyNode)
  self._vpcNodeMenu.insertItem(self.trUtf8('&Reaction'),self.vpcNewReactionNode)
  self._vpcNodeMenu.insertItem(self.trUtf8('&Filter'),self.vpcNewFilterNode)
  self._vpcNodeMenu.insertItem(self.trUtf8('Sin&k'),self.vpcNewSinkNode)

  populate(self._vpcView)

def populate(win):
  """ testing/demo code """
  from qtGui.Graph import Node

  n1 = Node.GraphNode(type='filter')
  filt = VPLib.VPFilterItem(win._canvas,graphParent=n1)
  filt.moveBy(50,200)
  filt.show()
  win.items.append(filt)

  n1 = Node.GraphNode(type='reaction')
  rxn = VPLib.VPReactionItem(win._canvas,
                            rxnFile='%s/qtGui/VisProg/exampleD/rxn1.cdxml'%RDConfig.RDCodeDir,
                            rxnImg='%s/qtGui/VisProg/exampleD/rxn1.gif'%RDConfig.RDCodeDir,
                             graphParent=n1)
  rxn.moveBy(200,125)
  rxn.show()
  win.items.append(rxn)

  n1 = Node.GraphNode(type='supply')
  supply = VPLib.VPSupplyItem(win._canvas,graphParent=n1)
  supply.moveBy(50,50)
  supply.show()
  win.items.append(supply)


  win.linkObjects(supply.getLinkPts()[-1],rxn.getLinkPts()[0])
  win.linkObjects(filt.getLinkPts()[-1],rxn.getLinkPts()[1])
  

