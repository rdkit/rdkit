# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time
from forms.itemdialog import ItemEdit as _Form

class ItemEdit(_Form):
  def __init__(self,name,
               parent=None,modal=0,fl=0):
    _Form.__init__(self,parent=parent,modal=modal,fl=fl)
    self.textBox.setText(name)
    self._taken = 1
    
  def updateFromItem(self,item):
    self.textBox.setText(item.attribs()['text'])

  def updateItem(self,item):
    if not self._taken:
      return
    lV = item.listView()
    lV.setNeedsSave(1)
    itemTxt = self.textBox.text()
    item.attribs()['text'] = itemTxt
    try:
      itemCol = lV.rawNameToSection('text')
    except ValueError:
      itemCol = -1
    else:
      item.setText(itemCol,itemTxt)

  def accept(self):
    self._taken = 1
    _Form.accept(self)

  def reject(self):
    self._taken = 0
    _Form.reject(self)
