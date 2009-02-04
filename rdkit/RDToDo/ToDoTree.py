# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time,weakref
from mx.DateTime import DateTimeFrom
from xml.dom import minidom
import OutlUtils,ToDoUtils
from ToDoTreeItem import ToDoTreeItem,ToDoTreeDupeItem
from TaskEditImpl import TaskEdit
from OutlineTree import OutlineTree
from ToDoListItem import ToDoListItem

class ToDoTree(OutlineTree):
  def __init__(self,parent,name,owner=None):
    OutlineTree.__init__(self,parent,name)
    self._childType = ToDoTreeItem
    self._owner = owner
    self._pixmapSize=20
    
  def _initColSettings(self):
    self._displayHeadings = {'text':'Item','due':'Due','priority':'Pri'}
    self._colOrder=['Pri','Due','Item',]
    self._canonColOrder=['Item','Pri','Due']
    self._revHeadings = {}
    for key in self._displayHeadings.keys():
      self._revHeadings[self._displayHeadings[key]] = key
    self._colWidths = [-1]*len(self._revHeadings)

  def _statusRefresh(self,start=None):
    if start is None:
      start = self
      node = start.firstChild()
      while node:
        node.calcStatus()
        node = node.nextSibling()
    else:
      try:
        start.setDirty(1)
        start.ultimateParent().calcStatus()
      except AttributeError:
        pass
      #start.calcStatus()
  def pixmapSize(self):
    return self._pixmapSize
  def setPixmapSize(self,size):
    self._pixmapSize=size

  def refresh(self,start=None):
    print 'refresh:',start
    self._statusRefresh(start=start)
    if start is None:
      self.updateVisibility()
      
  def newNode(self,name):
    node = OutlineTree.newNode(self,name)
    dateTxt = time.strftime('%d %B %Y',time.localtime())
    node.setAttribute('dateCreated',dateTxt)
    return node
  
  def addDomNode(self,node,where=None,after=None,select=0):
    if node.getAttribute('duplicate'):
      klass=ToDoTreeDupeItem
      dupe=True
      # make sure we aren't putting a duplicate under
      # its controller (that's very bad):
      if where:
        dupeId = node.getAttribute('duplicate')
        tmp=where
        while tmp:
          if tmp.GUID()==dupeId:
            print "cannot add a dupe under its controller"
            return
          tmp = tmp.parent()

    else:
      klass=self._childType
      dupe=False
    itm = OutlineTree.addDomNode(self,node,where=where,after=after,select=select,
                                 klass=klass)
    if where is not None:
      where.setDirty(1)
    
    percent = float(itm.attribs().get('percentDone',0.0))
    itm.setPixmap(0,ToDoUtils.GenCompletionPixmap(percent))
    #FIX: this is a STINKING hack
    if self._owner:
      #self._owner.listView.insertItem(itm)
      if percent != 100 and node.firstChild is None and not dupe:
        cpy = OutlineTree.addDomNode(self,node,klass=ToDoListItem)
        cpy.setPixmap(0,ToDoUtils.GenCompletionPixmap(percent))
        self.takeItem(cpy)
        self._owner.listView.insertItem(cpy)
        cpy._register()
        itm._partner = weakref.proxy(cpy)
        cpy._partner = weakref.proxy(itm)
    if dupe:
      itm.attribs()['duplicate']=node.getAttribute('duplicate')

    return itm

  def contextMenuEvent(self,evt):
    pos = QPoint(evt.pos().x(),evt.pos().y()-self.header().height())
    itm = self.itemAt(pos)
    sect = self.header().sectionAt(pos.x())
    if itm:
      itm._idMap = {}
      menu = None
      if not itm.attribs().get('duplicate',False):
        if sect == self.dispNameToSection('Pri'):
          menu = QPopupMenu(self)
          for i in range(10):
            id = menu.insertItem(str(i),itm._setPriFromPopup)
            itm._idMap[id] = i

        elif sect == self.dispNameToSection('Item'):
          menu = QPopupMenu(self)
          id = menu.insertItem(self.trUtf8('&Edit'),itm.editDialog)
          id = menu.insertSeparator()
          if itm.firstChild() is None:
            id = menu.insertItem(QIconSet(ToDoUtils.GenCompletionPixmap(100)),
                                 self.trUtf8('&Finish'),itm._updateStatusFromPopup)
            itm._idMap[id] = 100
            id = menu.insertSeparator()
          if itm.firstChild() is None:
            for i in range(11):
              if i != 10:
                lab = '&'
              else:
                lab = ''
              lab += '%d%%'%(i*10)
              percent = i*10
              id = menu.insertItem(QIconSet(ToDoUtils.GenCompletionPixmap(percent)),
                                   lab,
                                   itm._updateStatusFromPopup)
              itm._idMap[id] = percent

        elif sect == self.dispNameToSection('Due'):
          dueText = str(itm.text(sect))
          if len(dueText) > 0:
            menu = QPopupMenu(self)
            menu.insertItem(self.trUtf8('&Who set this date?'),itm._activateDateSetter)
      elif hasattr(itm,'getController'):
        if itm.getController():
          menu = QPopupMenu(self)
          menu.insertItem(self.trUtf8('&Activate Controller'),
                          itm._activateController)
          
      if menu:
        menu.exec_loop(self.mapToGlobal(pos))
        evt.consume()
        evt.accept()
      else:
        evt.ignore()
    else:
      evt.ignore()

  def editDialog(self,item):
    try:
      dlg = item._editDlg
    except:
      name = item.attribs()['text']
      try:
        priority = int(item.attribs()['priority'])
      except:
        priority = None
      dateStr = str(item.attribs().get('due',''))
      date = None
      if dateStr:
        try:
          date = DateTimeFrom(dateStr).tuple()
        except ValueError:
          pass
      try:
        note = item.attribs()['note']
      except:
        note = None
      dlg = TaskEdit(name,priority,date,note)
      item._editDlg = dlg
    else:
      dlg.updateFromItem(item)
    dlg.raiseW()
    res = dlg.exec_loop()
    if res:
      self.setNeedsSave(1)
      dlg.updateItem(item)
      if hasattr(item,'_partner'):
        dlg.updateItem(item._partner)
    self.refresh()

  def returnPress(self,item):
    OutlineTree.returnPress(self,item)
    self.refresh()

  def keyPressEvent(self,evt):
    handled = 0
    code = evt.key()
    item = None
    if evt.state() and Qt.ALT:
      item = self.selectedItem()
      if item and item.firstChild() is None:
        if code == Qt.Key_D:
          item._updateStatus(100)
          handled=1
        elif code == Qt.Key_0:
          item._updateStatus(0)
          handled=1
        elif code == Qt.Key_1:
          item._updateStatus(10)
          handled=1
        elif code == Qt.Key_2:
          item._updateStatus(20)
          handled=1
        elif code == Qt.Key_3:
          item._updateStatus(30)
          handled=1
        elif code == Qt.Key_4:
          item._updateStatus(40)
          handled=1
        elif code == Qt.Key_5:
          item._updateStatus(50)
          handled=1
        elif code == Qt.Key_6:
          item._updateStatus(60)
          handled=1
        elif code == Qt.Key_7:
          item._updateStatus(70)
          handled=1
        elif code == Qt.Key_8:
          item._updateStatus(80)
          handled=1
        elif code == Qt.Key_9:
          item._updateStatus(90)
          handled=1
    elif code==Qt.Key_P:
      self.verboseInfo()
      handled=1
    if not handled:
      OutlineTree.keyPressEvent(self,evt)
    else:
      if item:
        self.setNeedsSave(1)
        dateTxt = time.strftime('%d %B %Y',time.localtime())
        # FIX: implementation specific
        item._attribs['lastDoneModification']=dateTxt
        item.setDirty(1)
        if self._owner:
          if hasattr(item,'_partner'):
            item._partner._updateStatus(float(item.attribs()['percentDone']))
      self.refresh()
      evt.accept()


  def verboseInfo(self):
    node = self.firstChild()
    level = 0
    stack = [(level,node)]
    while len(stack)>0:
      level,node = stack.pop(0)
      print '  '*level,node.attribs().get('text'),node.GUID()
      child = node.firstChild()
      while child:
        stack.insert(0,(level+1,child))
        child = child.nextSibling()
      sib = node.nextSibling()
      while sib:
        stack.append((level,sib))
        sib = sib.nextSibling()


        
    
