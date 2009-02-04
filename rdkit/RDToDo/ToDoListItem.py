# $Id$
#
#  Copyright (C) 2003-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
from mx.DateTime import DateTimeFrom
from ToDoTreeItem import ToDoTreeItem
import time,weakref

class ToDoListItem(ToDoTreeItem):
  def __init__(self,*args,**kwargs):
    ToDoTreeItem.__init__(self,*args,**kwargs)
    self._partner = None
  def _register(self):
    QObject.connect(self.listView(),SIGNAL("itemRenamed(QListViewItem *,int)"),
                    self.itemRenamed)

  def partner(self):
    if self._partner is None:
      itm = self.listView()._tree.findItem('id',self.GUID())
      if itm:
        self._partner = weakref.proxy(itm)
      else:
        self._partner = itm
    return self._partner

  def itemRenamed(self,itm,col):
    if itm is not self: return
    if col == 0:
      txt = str(itm.text(col))
      itm._attribs['text'] = txt
      partner = self.partner()
      partner._attribs['text']=txt
      partner.setText(0,txt)
      lv = partner.listView()
      if lv is not None:
        lv.setNeedsSave(1)

  def key(self,col,ascending):
    if col==self.listView().dispNameToSection('Item'):
      return self._attribs['text']
    elif col==self.listView().dispNameToSection('Pri'):
      return str(self._attribs.get('priority',''))
    elif col==self.listView().dispNameToSection('Due'):
      d1 = self._attribs.get('due','')
      if d1:
        v1 = DateTimeFrom(str(d1)).ticks()
      else:
        v1 = 1e15
      return str(v1)
  def startRename(self,col):
    QListViewItem.startRename(self,col)

  def _updateStatusFromPopup(self,id):
    self.setDirty(1,doParents=0)
    percent = self._idMap[id]
    self._updateStatus(percent)

    # find the corresponding item in the tree:
    itm = self.partner()
    if itm:
      itm.setDirty(1,doParents=1)
      itm._updateStatus(percent)
      itm.listView().refresh()
      dateTxt = time.strftime('%d %B %Y',time.localtime())
      itm._attribs['lastDoneModification']=dateTxt

