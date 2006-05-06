# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ListSelectorDialog.ui'
#
# Created: Mon Jan 23 08:34:47 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class ListSelectorDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if not name:
      self.setName("ListSelectorDialog")


    ListSelectorDialogLayout = QVBoxLayout(self,4,2,"ListSelectorDialogLayout")

    self.listBox = QListBox(self,"listBox")
    self.listBox.setFrameShape(QListBox.StyledPanel)
    self.listBox.setFrameShadow(QListBox.Sunken)
    self.listBox.setSelectionMode(QListBox.Multi)
    ListSelectorDialogLayout.addWidget(self.listBox)

    layout5 = QHBoxLayout(None,0,2,"layout5")

    self.okButton = QPushButton(self,"okButton")
    layout5.addWidget(self.okButton)
    spacer1 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout5.addItem(spacer1)

    self.cancelButton = QPushButton(self,"cancelButton")
    layout5.addWidget(self.cancelButton)
    ListSelectorDialogLayout.addLayout(layout5)

    self.languageChange()

    self.resize(QSize(326,260).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
    self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)


  def languageChange(self):
    self.setCaption(self.__tr("Select"))
    self.okButton.setText(self.__tr("OK"))
    self.cancelButton.setText(self.__tr("Cancel"))


  def __tr(self,s,c = None):
    return qApp.translate("ListSelectorDialog",s,c)
