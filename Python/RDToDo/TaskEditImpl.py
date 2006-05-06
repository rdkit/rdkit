# $Id: TaskEditImpl.py 4899 2006-01-12 16:20:06Z glandrum $
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time
from mx.DateTime import DateTimeFrom
from forms.taskdialog import TaskEdit as _Form

class TaskEdit(_Form):
  def __init__(self,name,priority=None,date=None,note=None,
               parent=None,modal=0,fl=0):
    _Form.__init__(self,parent=parent,modal=modal,fl=fl)
    self.nameEdit.setText(name)
    if priority is not None and priority != '':
      self.priorityButton.setChecked(1)
      self.priorityBox.setEnabled(1)
      if priority < self.priorityBox.minValue():
        self.priorityBox.setMinValue(priority)
      if priority > self.priorityBox.maxValue():
        self.priorityBox.setMaxValue(priority)
      self.priorityBox.setValue(priority)

    if date is not None and date != '':
      self.dueButton.setChecked(1)
      self.dueBox.setEnabled(1)
    else:
      date = time.localtime()
    qD = QDate(date[0],date[1],date[2])
    self.dueBox.setDate(qD)

    if note is not None and note != '':
      self.note_edit.setText(note)
    self._taken = 1
    
  def updateFromItem(self,item):
    self.nameEdit.setText(item.attribs()['text'])

    try:
      pri = item.attribs()['priority']
    except KeyError:
      pass
    else:
      self.priorityBox.setValue(int(pri))
      
    dateStr = str(item.attribs().get('due',''))
    try:
      date = DateTimeFrom(dateStr).tuple()
    except ValueError:
      pass
    else:
      self.dueBox.setDate(QDate(date[0],date[1],date[2]))

  def updateItem(self,item):
    if not self._taken:
      return
    
    lV = item.listView()
    #lV.setNeedsSave(1)
    itemTxt = self.nameEdit.text()
    item.attribs()['text'] = itemTxt
    try:
      itemCol = lV.rawNameToSection('text')
    except ValueError:
      itemCol = -1
    else:
      item.setText(itemCol,itemTxt)
    if self.dueButton.isChecked():
      qD = self.dueBox.date()
      year = qD.year()
      month = qD.month()
      day = qD.day()
      l = [0]*9
      l[0] = year
      l[1] = month
      l[2] = day
      dateTxt = time.strftime('%d %B %Y',tuple(l))
    else:
      dateTxt = ''
    item.attribs()['due'] = dateTxt
    try:
      dueCol = lV.rawNameToSection('due')
    except ValueError:
      dueCol = -1
    else:  
      item.setText(dueCol,dateTxt)

    if self.priorityButton.isChecked():
      pri = self.priorityBox.value()
      item.attribs()['priority'] = pri
      try:
        priCol = lV.rawNameToSection('priority')
      except ValueError:
        priCol = -1
      else:  
        item.setText(priCol,str(pri))
    txt = self.note_edit.text()
    item.attribs()['note'] = txt
    item.setDirty(1)
    p = item.parent()
    while p:
      p.setDirty(1)
      p = p.parent()
  #
  #
  #

  def dueClick(self):
    if self.dueButton.isChecked():
      self.dueBox.setEnabled(1)
    else:
      self.dueBox.setEnabled(0)
      
  def priorityClick(self):
    if self.priorityButton.isChecked():
      self.priorityBox.setEnabled(1)
    else:
      self.priorityBox.setEnabled(0)

  def accept(self):
    self._taken = 1
    _Form.accept(self)

  def reject(self):
    self._taken = 0
    _Form.reject(self)

