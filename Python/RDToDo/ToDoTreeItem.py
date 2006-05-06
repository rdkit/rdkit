# $Id: ToDoTreeItem.py 4900 2006-01-12 18:06:52Z glandrum $
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time,weakref
from mx.DateTime import DateTimeFrom
import ToDoUtils
from icons import icons

from OutlineTreeItem import OutlineTreeItem

class ToDoTreeItem(OutlineTreeItem):
  """ an item for a ToDo tree

  attributes needed:
    ['text','priority','due','note','percentDone']

  """
  def __init__(self,*args,**kwargs):
    apply(OutlineTreeItem.__init__,(self,)+args,kwargs)
    self._dirty=1
    self._overDue=0
    self._duplicates={}

  def _setPriFromPopup(self,id):
    pri = self._idMap[id]
    self._attribs['priority']=str(pri)
    try:
      priCol = self.listView().rawNameToSection('priority')
    except ValueError:
      priCol = -1
    else:  
      self.setText(priCol,str(pri))
      
  def _updateStatus(self,percent):
    #print '  '*self.getLevel(),'u:',self._attribs['text'],
    #if hasattr(self,'_controller'):
    #  print '-dupe',
    #print
    self._attribs['percentDone']=str(percent)
    self.setPixmap(0,ToDoUtils.GenCompletionPixmap(percent))

    for dupe in self._duplicates.values():
      dupe()._updateStatus(percent)

  def _updateDate(self,date):
    try:
      col = self.listView().rawNameToSection('due')
    except ValueError:
      col = -1
    if date is None:
      if col != -1:
        pm = QPixmap()
        pm.loadFromData(icons.blankPng)
        self.setPixmap(col,pm)
      return
    #date = tuple([int(x) for x in date])
    if date and date[0] != 4223 and self.getPercentDone() != 100 and \
           time.mktime(date) - time.time() < 0:
      self._overDue = 1
    else:
      self._overDue = 0
    if col != -1:
      dateTxt = ToDoUtils.TextifyDate(date)
      if len(dateTxt)>0 and date != self.getDate():
        dateTxt = '(%s)'%(dateTxt)
      self.setText(col,dateTxt)
      pm = QPixmap()
      if self._overDue:
        pm.loadFromData(icons.exclaimPng)
      else:
        pm.loadFromData(icons.blankPng)
      self.setPixmap(col,pm)

  def _activateDateSetter(self):
    setter = self
    date = self.getDate()
    child = self.mostPressingChild()
    if child is not None:
      cDate = child.getDate()
      if date is None or (cDate is not None and cDate < date):
        setter = child
    p = setter.parent()
    while p:
      p.setOpen(1)
      p = p.parent()
    self.listView().setSelected(setter,1)  

  def getPercentDone(self):
    return int(self._attribs.get('percentDone'))
  
  def getLevel(self):
    level = 0
    p = self.parent()
    while p is not None:
      p = p.parent()
      level += 1
    return level  
      
  def getDate(self):
    """ returns this object's date as a python tuple,
    None if it's not set """
    d = self._attribs.get('due',None) 
    if d is None or len(d)==0:
      date = None
    else:
      dateTxt = str(d)
      date = tuple([int(x) for x in DateTimeFrom(dateTxt).tuple()])
    return date

  def mostPressingChild(self):
    #print '%smpc: '%('  '*self.getLevel()),self._attribs['text']
    # if we are done, we have no pressing children
    if self.getPercentDone() == 100:
      return None


    child = self.firstChild()
    res = child
    d1 = None
    while child:
      # no point considering finished children
      if child.getPercentDone() != 100:
        d = child.getDate()
        #print '%s->'%('  '*self.getLevel()),child._attribs['text'],d
        if d is not None:
          if d1 is None or d < d1:
            d1 = d
            res = child
        pressing = child.mostPressingChild()
        if pressing is not None:
          d = pressing.getDate()
          if d is not None and (d1 is None or d < d1):
            d1 = d
            res = pressing
      child = child.nextSibling()  
    return res  
    
  
  def dirty(self):
    return self._dirty
  def setDirty(self,dirty=1,doParents=1):
    self._dirty=dirty
    p = self.parent()
    if doParents and p is not None:
      p.setDirty(1,doParents=doParents)
    for dupe in self._duplicates.values():
      dupe().setDirty(dirty=dirty,doParents=doParents)
  
  def _updateStatusFromPopup(self,id):
    #print '***************************************'
    self.setDirty(1,doParents=1)
    percent = self._idMap[id]
    self._updateStatus(percent)
    self.listView().refresh()

    dateTxt = time.strftime('%d %B %Y',time.localtime())
    self._attribs['lastDoneModification']=dateTxt

    if hasattr(self,'_partner'):
      self._partner._updateStatus(percent)
    
  def calcStatus(self,useChildren=1):
    """

      returns a (percent,date) tuple

      NOTE:  calcStatus
        *updates* the internal done percent and
        *calculates* the due date (the display only is updated)


    """
    #print '  '*self.getLevel(),self._attribs['text'],
    #if hasattr(self,'_controller'):
    #  print '-dupe',
    #print self.dirty()
    if not self.dirty():
      percent=int(self._attribs.get('percentDone',0))
      date = self.getDate()
      #print 'undirty date:',self._attribs['text'],repr(date)
      if date is None:
        child = self.mostPressingChild()
        if child:
          date = child.getDate()
        else:
          date = None
      if date:
        self._updateDate(date)
      return percent,date      

    child = self.firstChild()
    # this seems a bit hacky, but we know how python
    #  compares tuples, so it'll be ok
    earlyDate = self.getDate()
    if earlyDate is None:
      earlyDate = (4223,0,0,0,0,0,0,0,0)

    if useChildren and child is not None:
      nChildren = 1
      accum,date = child.calcStatus(useChildren=1)
      if date is not None and accum < 100 and date < earlyDate:
        earlyDate = date
      else:
        pass
      child = child.nextSibling()
      while child:
        nChildren += 1
        percent,date = child.calcStatus(useChildren=1)
        if date is not None and percent < 100 and date < earlyDate:
          earlyDate = date
        accum += percent
        child = child.nextSibling()
    else:
      # we're okay, we don't need to do anything
      if hasattr(self,'_controller'):
        accum,earlyDate = self.getController().calcStatus()
      else:
        accum = int(self._attribs.get('percentDone',0))
        earlyDate = self.getDate()
      nChildren = 1

    percent = accum/nChildren
    self.setDirty(0,doParents=0)
    self._updateStatus(percent)
    self._updateDate(earlyDate)
    return percent,earlyDate

  def paintCell(self,painter,cg,col,width,align):
    if self._attribs.get('duplicate',0):
      needRestore=True
      painter.save()
      fnt = painter.font()
      fnt.setItalic(True)
      painter.setFont(fnt)
    else:
      needRestore=False
    QListViewItem.paintCell(self,painter,cg,col,width,align)

    if needRestore:
      painter.restore()

class ToDoTreeDupeItem(ToDoTreeItem):
  """ an item for a ToDo tree

  attributes needed:
    ['duplicate']

  """
  def __init__(self,*args,**kwargs):
    apply(ToDoTreeItem.__init__,(self,)+args,kwargs)
    self._controller=None
    for i in range(self.listView().columns()):
      self.setRenameEnabled(i,False)
    self.setDropEnabled(False)

  def acceptsChildren(self):
    return False

  def getController(self):
    if not self._controller:
      node = self._findController()
      if node:
        self._controller = weakref.ref(node)
        node._duplicates[self.GUID()]=weakref.ref(self)
    return self._controller()
  def _findController(self,toDo=None):
    dupId = self._attribs['duplicate']

    if toDo is None:
      toDo = self.listView().firstChild()
    while toDo:
      if toDo._attribs.get('id','')==dupId:
        return toDo
      else:
        child = toDo.firstChild()
        while child:
          r = self._findController(toDo=child)
          if r:
            return r
          child = child.nextSibling()
      toDo = toDo.nextSibling()
    return None
    
  def getPercentDone(self):
    ctrl = self.getController()
    assert ctrl,'No Controller found'
    if ctrl:
      return ctrl.getPercentDone()
    
  def getDate(self):
    ctrl = self.getController()
    assert ctrl,'No Controller found'
    if ctrl:
      return ctrl.getDate()

  def activate(self):
    return
  def _activateController(self):
    ctrl = self.getController()
    if not ctrl:
      return
    p = ctrl.parent()
    while p:
      p.setOpen(True)
      p = p.parent()
    self.listView().setSelected(ctrl,True)
    return False

