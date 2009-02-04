# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time
import ToDoUtils
from utils import GUIDs


class OutlineTreeItem(QListViewItem):
  """ a generic item to go in an outline tree


   attributes used:
      ['text']

  """
  def __init__(self,*args,**kwargs):
    apply(QListViewItem.__init__,(self,)+args)
    self._domNode = None
    self._attribs = {}
    self._maxLineLen=80
    self._attach = None
    self.setRenameEnabled(0,1)
    parent = args[0]
    QObject.connect(self.listView(),SIGNAL("itemRenamed(QListViewItem *,int)"),
                    self.itemRenamed)
    self._idMap = {}

  def __str__(self):
    return self.toText()

  def toText(self,depth=0):
    res = '  '*depth+str(self._attribs.get('text'))+'\n'
    child = self.firstChild()
    while child is not None:
      res += child.toText(depth=depth+1)
      child = child.nextSibling()
    return res

  def node(self):
    return self._domNode
  def setNode(self,node):
    self._domNode = node

  def attribs(self):
    return self._attribs
  def clearAttribs(self):
    self._attribs = {}

  def collapseChildren(self,state=0,recurse=1):
    child = self.firstChild()
    while child:
      if recurse:
        child.collapseChildren(state=state,recurse=1)
      child = child.nextSibling()  
      self.setOpen(state)
        
  def setText(self,col,itemText):
    itemText = str(itemText)
    if itemText.find('\n') > 0:
      itemText = itemText[:itemText.find('\n')]
    #if len(itemText) > self._maxLineLen:
    #  itemText = itemText[:self._maxLineLen]
    QListViewItem.setText(self,col,itemText)
      
  def setGUID(self,guid=None):
    if guid is None:
      guid = GUIDs.getGUID()
    self._attribs['id']=guid
  def GUID(self):
    return self._attribs.get('id','')
  def ultimateParent(self):
    """ returns the top-most parent (or ourself) """
    p = self
    while p.parent():
      p = p.parent()
    return p  
  def acceptsChildren(self):
    return True

  #
  #  signals and slots and handlers and stuff
  #
  def itemRenamed(self,itm,col):
    if itm is not self: return
    if col == 0:
      txt = str(itm.text(col))
      itm._attribs['text'] = txt
    lv = self.listView()
    if lv is not None:
      lv.setNeedsSave(1)
    if hasattr(self,'_partner'):
      self._partner.setText(col,txt)
  def editDialog(self):
    self.listView().editDialog(self)
  def startRename(self,col):
    self.listView().setNeedsSave(1)
    QListViewItem.startRename(self,col)

