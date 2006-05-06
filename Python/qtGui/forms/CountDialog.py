# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'CountDialog.ui'
#
# Created: Mon Jan 23 08:34:47 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class CountDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if not name:
      self.setName("CountDialog")


    CountDialogLayout = QVBoxLayout(self,2,4,"CountDialogLayout")

    self.labelText = QLabel(self,"labelText")
    CountDialogLayout.addWidget(self.labelText)

    layout11 = QHBoxLayout(None,0,6,"layout11")

    self.boxLabel = QLabel(self,"boxLabel")
    self.boxLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    layout11.addWidget(self.boxLabel)

    self.countBox = QSpinBox(self,"countBox")
    self.countBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.countBox.sizePolicy().hasHeightForWidth()))
    self.countBox.setMinimumSize(QSize(75,0))
    layout11.addWidget(self.countBox)
    spacer6 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout11.addItem(spacer6)
    CountDialogLayout.addLayout(layout11)

    layout9 = QHBoxLayout(None,0,6,"layout9")

    self.okButton = QPushButton(self,"okButton")
    self.okButton.setDefault(1)
    layout9.addWidget(self.okButton)
    spacer5 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout9.addItem(spacer5)

    self.cancelButton = QPushButton(self,"cancelButton")
    layout9.addWidget(self.cancelButton)
    CountDialogLayout.addLayout(layout9)

    self.languageChange()

    self.resize(QSize(354,118).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
    self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)


  def languageChange(self):
    self.setCaption(self.__tr("Count Dialog"))
    self.labelText.setText(self.__tr("textLabel2"))
    self.boxLabel.setText(self.__tr("Count:"))
    self.okButton.setText(self.__tr("OK"))
    self.cancelButton.setText(self.__tr("Cancel"))


  def __tr(self,s,c = None):
    return qApp.translate("CountDialog",s,c)
