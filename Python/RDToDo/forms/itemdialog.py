# Form implementation generated from reading ui file 'itemdialog.ui'
#
# Created: Fri Feb 21 07:27:11 2003
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class ItemEdit(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("ItemEdit")

    self.resize(511,273)
    self.setCaption(self.trUtf8("Edit Item"))
    self.setSizeGripEnabled(1)

    ItemEditLayout = QVBoxLayout(self,4,4,"ItemEditLayout")

    self.textBox = QTextEdit(self,"textBox")
    QToolTip.add(self.textBox,self.trUtf8("enter the text for the item (only the first line is displayed)"))
    ItemEditLayout.addWidget(self.textBox)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setEnabled(0)
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer)

    self.buttonOk = QPushButton(self,"buttonOk")
    self.buttonOk.setText(self.trUtf8("OK"))
    self.buttonOk.setAccel(0)
    self.buttonOk.setAutoDefault(1)
    self.buttonOk.setDefault(1)
    Layout1.addWidget(self.buttonOk)

    self.buttonCancel = QPushButton(self,"buttonCancel")
    self.buttonCancel.setText(self.trUtf8("Cancel"))
    self.buttonCancel.setAccel(0)
    self.buttonCancel.setAutoDefault(1)
    Layout1.addWidget(self.buttonCancel)
    ItemEditLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))

  def dueClick(self):
    print "ItemEdit.dueClick(): Not implemented yet"

  def priorityClick(self):
    print "ItemEdit.priorityClick(): Not implemented yet"
