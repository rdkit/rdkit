# $Id: ToDoListView.py 4899 2006-01-12 16:20:06Z glandrum $
#
#  Copyright (C) 2003-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
from mx.DateTime import DateTimeFrom
import weakref,time
import ToDoUtils
from ToDoTree import ToDoTree
from TaskEditImpl import TaskEdit

class ToDoListView(QListView):
  def __init__(self,tree,*args,**kwargs):
    QListView.__init__(self,*args,**kwargs)
    self._initColSettings()
    self._tree = tree
               
  def _initColSettings(self):
    self._displayHeadings = {'text':'Item','due':'Due','priority':'Pri'}
    self._colOrder=['Pri','Due','Item',]
    self._canonColOrder=['Item','Pri','Due']
    self._revHeadings = {}
    for key in self._displayHeadings.keys():
      self._revHeadings[self._displayHeadings[key]] = key
    self._colWidths = [-1]*len(self._revHeadings)

  def dispNameToSection(self,name):
    return self._canonColOrder.index(name)
  
  def rawNameToSection(self,name):
    return self._canonColOrder.index(self._displayHeadings[name])
  
  def contextMenuEvent(self,evt):
    pos = QPoint(evt.pos().x(),evt.pos().y()-self.header().height())
    itm = self.itemAt(pos)
    sect = self.header().sectionAt(pos.x())
    if itm:
      itm._idMap = {}
      menu = None
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
      dlg.updateItem(item)
      partner = item.partner()
      print 'update:',partner
      if partner:
        dlg.updateItem(partner)
        partner.listView().refresh()


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
    elif code == Qt.Key_Space:
      item = self.selectedItem()
      if item is not None:
        item.startRename(0)
        handled = 1
          
    if not handled:
      QListView.keyPressEvent(self,evt)
    else:
      if item:
        self._tree.setNeedsSave(1)
        dateTxt = time.strftime('%d %B %Y',time.localtime())
        # FIX: implementation specific
        item._attribs['lastDoneModification']=dateTxt
        item.setDirty(1)
        if self._tree:
          partner = item.partner()
          if partner:
            partner._updateStatus(int(item.attribs()['percentDone']))
            partner.setDirty(1)
      self._tree.refresh()
      evt.accept()

