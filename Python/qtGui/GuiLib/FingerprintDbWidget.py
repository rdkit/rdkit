#
#  Copyright (C) 2003  Rational Discovery LLC
#    All Rights Reserved
#
"""  Defines a widget for working with fingerprints in databases

"""
from qt import *
from qtGui.DbWidget import insertDbWidget
from qtGui.DbQueryWidgetImpl import DbQueryWidget

def insertFingerprintDbWidget(parent,*args,**kwargs):
  widg = insertDbWidget(parent,FingerprintDbWidget,*args,**kwargs)
  return widg

class FingerprintDbWidget(DbQueryWidget):
  def __init__(self,parent,*args,**kwargs):
    DbQueryWidget.__init__(self,parent,*args,**kwargs)
    self._initFpBox()

  def _initFpBox(self):
    parent = self.parent

    # and here's the a widget for specifying fingerprint column:
    spacer_2 = QSpacerItem(10,60,QSizePolicy.Minimum,QSizePolicy.Expanding)
    parent.layout().addItem(spacer_2)

    layout = QHBoxLayout(None,2,2,'fpColLayout')
    self.fpColLabel = QLabel(parent,"fpColLabel")
    self.fpColLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    self.fpColLabel.setText(self.trUtf8('Fingerprint Column'))
    layout.addWidget(self.fpColLabel)
    QToolTip.add(self.fpColLabel,self.trUtf8("Specify which column in the results contains the fingerprints.  This can be left blank, but the resulting query will be slow."))

    self.fpColBox = QLineEdit(parent,"fpColBox")
    QToolTip.add(self.fpColBox,self.trUtf8("Specify which column in the results contains the fingerprints.  This can be left blank, but the resulting query will be slow."))
    layout.addWidget(self.fpColBox)

    self.fpColButton = QToolButton(parent,"fpColButton")
    self.fpColButton.setText(self.trUtf8("Choose"))
    QToolTip.add(self.fpColButton,self.trUtf8("Select the fingerprint column from a list"))
    self.connect(self.fpColButton,SIGNAL("clicked()"),self.fpColChooseClick)
    layout.addWidget(self.fpColButton)

    self.fpColLabel.setEnabled(0)
    self.fpColBox.setEnabled(0)
    self.fpColButton.setEnabled(0)

    self.insertFileClickCallback(self.enableFpColSelection)
    parent.layout().addItem(layout)

  def getFpColumn(self):
    return str(self.fpColBox.text())
  
  def enableFpColSelection(self):
    self.fpColLabel.setEnabled(1)
    self.fpColBox.setEnabled(1)
    self.fpColButton.setEnabled(1)
    
  def fpColChooseClick(self):
    initVal = str(self.getFpColumn())
    nameL = self.pickColumnNames(startingSelection=initVal,
                                 allowMultiple=0)
    if nameL:
      self.fpColBox.setText(str(nameL[0]))
