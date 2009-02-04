# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time
from xml.dom import minidom
import OutlUtils,ToDoUtils
from OutlineTreeItem import OutlineTreeItem
from ItemEditImpl import ItemEdit

class OutlineTree(QListView):
  def __init__(self,parent,name):
    QListView.__init__(self,parent,name)
    self.setSorting(-1)
    self.setAcceptDrops(1)
    self.setDefaultRenameAction(QListView.Accept)
    self._dragging = 0
    self._dragItem = None
    self._details={}
    self._childType = OutlineTreeItem
    self._needsSave = 0
    self._initColSettings()
    QObject.connect(self,SIGNAL("returnPressed(QListViewItem *)"),
                    self.returnPress)
    self._invisibleThings = {}
    
  def _initColSettings(self):
    self._displayHeadings = {'text':'Text'}
    self._colOrder=['Text']
    self._canonColOrder=['Text']
    self._revHeadings = {}
    for key in self._displayHeadings.keys():
      self._revHeadings[self._displayHeadings[key]] = key
    self._colWidths = [-1]*len(self._revHeadings)

  def _colInit(self):
    while self.columns():
      self.removeColumn(0)
    # place the columns, then move them about (don't get me started)
    for i in range(len(self._colOrder)):
      self.addColumn(self._canonColOrder[i],-1)

    hdr = self.header()
    for i in range(len(self._colOrder)):
      startIdx = self._canonColOrder.index(self._colOrder[i])
      hdr.moveSection(startIdx,i)
      hdr.resizeSection(startIdx,int(self._colWidths[i]))

  def updateVisibility(self,start=None):
    if start is None:
      start = self
    node = start.firstChild()
    while node:
      self.updateVisibility(start=node)
      nodeAttribs = node.attribs()
      visible = 1
      for attrib in self._invisibleThings.keys():
        try:
          if nodeAttribs[attrib] in self._invisibleThings[attrib]:
            visible=0
        except:
          pass
      node.setVisible(visible)  
      node = node.nextSibling()

  def dispNameToSection(self,name):
    return self._canonColOrder.index(name)
  
  def rawNameToSection(self,name):
    return self._canonColOrder.index(self._displayHeadings[name])
  
  def empty(self):
    node = self.firstChild()
    while node:
      sib = node.nextSibling()
      self.takeItem(node)
      node = sib

  def headings(self):
    hdr = self.header()
    res =['']*hdr.count()
    for i in range(hdr.count()):
      res[i]  = str(hdr.label(i))
    return res

  def collapseAll(self,state=0):
    child = self.firstChild()
    while child:
      child.collapseChildren(state=state,recurse=1)
      child.setOpen(state)
      child = child.nextSibling()
      
  def refresh(self,start=None):
    if start is None:
      self.updateVisibility(start=None)

  def newNode(self,name):
    doc = minidom.Document()
    node = doc.createElement('outline')
    node.setAttribute('text','-- Unknown --')
    return node

  def addDomNode(self,node,where=None,after=None,select=0,klass=None):
    if klass is None:
      klass = self._childType
    heads = self.headings()
    things = [None]*len(heads)
    for i in range(len(heads)):
      attr = self._revHeadings[heads[i]]
      things[i] = node.getAttribute(attr)
    if where is None:
      thing = self
    else:
      thing = where
    itm = klass(thing,*things)


    #itm = apply(self._childType,(thing,)+tuple(things))  
    itm.setDragEnabled(True)
    if itm.acceptsChildren():
      itm.setDropEnabled(True)
    else:
      itm.setDropEnabled(False)
    itm.setNode(node)
    itm.clearAttribs()
    attrs = node.attributes
    for i in range(attrs.length):
      nm = attrs.item(i).name
      val = attrs.item(i).value
      itm.attribs()[nm] = val
    if not itm.attribs().get('id',''):
      itm.setGUID()
    if select:
      self.setSelected(itm,1)

    node = node.firstChild
    lastNode = None
    while node:
      if node.nodeType == node.ELEMENT_NODE and node.nodeName == 'outline':
        lastNode = self.addDomNode(node,where=itm,after=lastNode)
      node = node.nextSibling
    if after:
      itm.moveItem(after)
    return itm
  
  def addOPML(self,data,where=None,new=0):
    self.setNeedsSave(1)
    dom = minidom.parseString(data)
    if new:
      self._details = {}
      self._initColSettings()
      head = dom.getElementsByTagName('head')[0]
      for node in head.childNodes:
        nm = node.nodeName
        if node.firstChild:
          val = node.firstChild.data
          if nm == 'colOrder':
            self._colOrder = val.split(',')
            try:
              while 1:
                self._colOrder.remove('')
            except ValueError:
              pass
          elif nm == 'colWidths':
            self._colWidths = val.split(',')
            try:
              while 1:
                self._colWidths.remove('')
            except ValueError:
              pass
          else:
            self._details[nm] = val
      self._colInit()
    root = dom.getElementsByTagName('body')[0]
    node = root.firstChild
    lastNode = None
    while node:
      if node.nodeType == node.ELEMENT_NODE and node.nodeName == 'outline':
        lastNode=self.addDomNode(node,where=where,after=lastNode)
      node = node.nextSibling
    self.refresh(start=where)
    return lastNode
  
  def _delItem(self,item):
    parent = item.parent()
    if parent is None:
      parent = self
    # FIX: hack-hack-hack
    if hasattr(item,'_duplicates'):
      for dupe in item._duplicates.values():
        dupe = dupe()
        dupe._controller=None
        if hasattr(dupe,'setDirty'):
          dupe.setDirty(True,doParents=True)
    item.setGUID('')
    parent.takeItem(item)
    self.refresh(start=parent)

  def toxml(self,item=None,clearIDs=0):
    doc = minidom.Document()
    opml = doc.createElement('opml')
    opml.setAttribute('version','1.0')
    head = doc.createElement('head')
    self._details['dateModified']='%s GMT'%time.asctime(time.gmtime())
    for key in self._details.keys():
      n = doc.createElement(key)
      txt = doc.createTextNode(str(self._details[key]))
      n.appendChild(txt)
      head.appendChild(n)

    hdr = self.header()
    colNameStr = ''
    colWidthStr = ''
    for i in range(hdr.count()):
      section = hdr.mapToSection(i)
      colNameStr = '%s%s,'%(colNameStr,str(hdr.label(section)))
      colWidthStr = '%s%s,'%(colWidthStr,str(hdr.sectionSize(section)))
    colNameStr = colNameStr[:-1]
    n = doc.createElement('colOrder')
    txt = doc.createTextNode(colNameStr)
    n.appendChild(txt)
    head.appendChild(n)

    colWidthStr = colWidthStr[:-1]
    n = doc.createElement('colWidths')
    txt = doc.createTextNode(colWidthStr)
    n.appendChild(txt)
    head.appendChild(n)
    
    opml.appendChild(head)
    body = doc.createElement('body')
    if item is not None:
      node = OutlUtils.outlItemToDOM(item,doc=doc,parentNode=body,clearIDs=clearIDs)
      body.appendChild(node)
    else:
      item = self.firstChild()
      while item:
        node = OutlUtils.outlItemToDOM(item,doc=doc,parentNode=body,clearIDs=clearIDs)
        body.appendChild(node)
        item = item.nextSibling()
        
    opml.appendChild(body)
    doc.appendChild(opml)
    return doc.toxml()
    
  def findSiblingOver(self,itm):
    if itm is None:
      return None
    parent = itm.parent()
    if parent is None:
      parent = itm.listView()
    res = None
    child = parent.firstChild()
    while child and child != itm:
      sib = child.nextSibling()
      if sib == itm:
        res = child
      child = sib
    return res

  def needsSave(self):
    return self._needsSave
  def setNeedsSave(self,needs=1):
    self._needsSave = needs

  def findItem(self,field,val,substring=0,caseSensitive=1):
    consider = [self.firstChild()]
    res = None
    while len(consider):
      entry = consider.pop()
      if OutlUtils.findHelper(entry,field,val,substring=substring,
                              caseSensitive=caseSensitive):
        res = entry
        break
      else:
        next = entry.firstChild()
        while next:
          consider.append(next)
          next = next.nextSibling()
    return res      
          
    
  #    
  #
  #
  def startDrag(self):
    item = self.selectedItem()
    self._dragItem = item
    dragTxt = self.toxml(item=item)
    d = QTextDrag(dragTxt,self)
    d.dragCopy()
    
  def contextMenuEvent(self,evt):
    pos = QPoint(evt.pos().x(),evt.pos().y()-self.header().height())
    itm = self.itemAt(pos)
    sect = self.header().sectionAt(pos.x())
    if itm:
      menu = QPopupMenu(self)
      id = menu.insertItem(self.trUtf8('&Edit'),itm.editDialog)

      if menu:
        res = menu.exec_loop(self.mapToGlobal(pos))
        evt.consume()
        evt.accept()
      else:
        evt.ignore()
    else:
      evt.ignore()
    
  def dragEnterEvent(self,evt):
    evt.accept(QTextDrag.canDecode(evt))
    self._dragging = 1

  def dropEvent(self,evt):
    # FIX: we should be able to drop raw text onto an item
    txt = QString()
    res = QTextDrag.decode(evt,txt)
    txt = str(txt)
    if res:
      pos = QPoint(evt.pos().x(),evt.pos().y()-self.header().height())
      item = self.itemAt(pos)
      
      ok = True
      self._dragging = 0
      if item:
        if item.acceptsChildren():
          parent = item
          while parent:
            if parent == self._dragItem:
              ok = False
              break
            parent = parent.parent()
        else:
          ok = False
      else:
        item = None
      if ok:
        self.setNeedsSave(1)
        added=self.addOPML(txt,item)
        if added and self._dragItem:
          # FIX: hack-hack-hack
          if hasattr(self._dragItem,'_duplicates'):
            dupes = self._dragItem._duplicates.values()
          else:
            dupes = []
          self._delItem(self._dragItem)
          for dupe in dupes:
            dupe().getController()
        if item is not None:
          item.setOpen(1)

  def editDialog(self,item):
    print 'outleditdialog'
    try:
      dlg = item._editDlg
    except:
      text = item.attribs()['text']
      dlg = ItemEdit(text)
      item._editDlg = dlg
    else:
      dlg.updateFromItem(item)
    dlg.raiseW()
    dlg.exec_loop()
    dlg.updateItem(item)
    self.refresh()

  def keyPressEvent(self,evt):
    handled = 0
    code = evt.key()
    # the switching here is kind of putrid because of some windows
    #  weirdness I need to track down
    item = None
    needRefresh = []
    if code == Qt.Key_Tab and not(evt.state() and Qt.SHIFT):
      item = self.selectedItem()
      if item is not None:
        parent = item.parent()
        # move down the hierarchy
        if parent is None:
          parent = item.listView()
        overItm = self.findSiblingOver(item)
        if overItm and overItm.acceptsChildren():
          parent.takeItem(item)
          needRefresh.append(parent)
          lastChild = None
          if overItm.childCount() > 0:
            c = overItm.firstChild()
            while c.nextSibling():
              c = c.nextSibling()
            lastChild = c
          overItm.insertItem(item)
          if lastChild:
            item.moveItem(lastChild)
          overItm.setOpen(1)
          needRefresh.append(overItm)
          handled = 1
        item.listView().setSelected(item,1)
    elif code in [Qt.Key_BackTab,Qt.Key_Tab]:
      item = self.selectedItem()
      if item is not None:
        parent = item.parent()
        # move back up the hierarchy
        if parent is not None:
          # our parent isn't the top, so let's move up one
          overItm = parent.parent()
          if overItm is None:
            overItm = item.listView()
            needRefresh.append(parent)
          else:
            needRefresh.append(parent)
            needRefresh.append(overItm)
          parent.takeItem(item)
          overItm.insertItem(item)
          item.moveItem(parent)
          try:
            overItm.setOpen(1)
          except TypeError:
            pass
          handled = 1
        item.listView().setSelected(item,1)
    elif code == Qt.Key_Up and (evt.state() and Qt.SHIFT):
      item = self.selectedItem()
      if item is not None:
        overItm = self.findSiblingOver(item)
        # move back up the hierarchy
        if overItm is not None:
          overoverItm = self.findSiblingOver(overItm)
          if overoverItm:
            item.moveItem(overoverItm)
          else:
            if item.parent():
              parent = item.parent()
            else:
              parent = item.listView()
            parent.takeItem(item)
            parent.insertItem(item)
        handled = 1
        item.listView().setSelected(item,1)
    elif code == Qt.Key_Down and (evt.state() and Qt.SHIFT):
      item = self.selectedItem()
      if item is not None:
        underItm = item.nextSibling()
        if underItm is not None:
          item.moveItem(underItm)
        handled = 1
        item.listView().setSelected(item,1)
    elif code == Qt.Key_Space:
      item = self.selectedItem()
      if item is not None:
        item.startRename(0)
        handled = 1
      
    if not handled:
      QListView.keyPressEvent(self,evt)
    else:
      if len(needRefresh):
        self.setNeedsSave(1)
      for entry in needRefresh:
        try:
          entry.setDirty(1)
        except AttributeError:
          pass
      self.refresh()
      evt.accept()


  def returnPress(self,item):
    self.setNeedsSave(1)
    node = self.newNode('-- Unnamed --')
    newItm = self.addDomNode(node,where=item.parent(),after=item,
                             select=1)
    newItm.startRename(0)
    self.ensureItemVisible(newItm)
      
      

