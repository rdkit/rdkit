# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AtomCorrespondenceDialog.ui'
#
# Created: Thu Oct 19 08:53:57 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.16
#
# WARNING! All changes made in this file will be lost!


from qt import *


class AtomCorrespondenceDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if not name:
      self.setName("AtomCorrespondenceDialog")

    self.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.sizePolicy().hasHeightForWidth()))
    self.setMinimumSize(QSize(0,0))

    AtomCorrespondenceDialogLayout = QVBoxLayout(self,11,6,"AtomCorrespondenceDialogLayout")

    layout4 = QHBoxLayout(None,0,6,"layout4")

    self.showCorrespondenceBox = QCheckBox(self,"showCorrespondenceBox")
    self.showCorrespondenceBox.setChecked(1)
    layout4.addWidget(self.showCorrespondenceBox)

    self.makeCorrespondenceButton = QPushButton(self,"makeCorrespondenceButton")
    layout4.addWidget(self.makeCorrespondenceButton)

    self.clearCorrespondenceButton = QPushButton(self,"clearCorrespondenceButton")
    layout4.addWidget(self.clearCorrespondenceButton)
    AtomCorrespondenceDialogLayout.addLayout(layout4)

    self.refineButton = QPushButton(self,"refineButton")
    AtomCorrespondenceDialogLayout.addWidget(self.refineButton)

    layout3 = QGridLayout(None,1,1,0,6,"layout3")

    self.textLabel3 = QLabel(self,"textLabel3")
    self.textLabel3.setAlignment(QLabel.AlignCenter)

    layout3.addWidget(self.textLabel3,1,2)

    self.rmsdBox = QLineEdit(self,"rmsdBox")
    self.rmsdBox.setReadOnly(1)

    layout3.addWidget(self.rmsdBox,0,2)

    self.shapeScoreBox = QLineEdit(self,"shapeScoreBox")
    self.shapeScoreBox.setReadOnly(1)

    layout3.addWidget(self.shapeScoreBox,0,1)

    self.energyBox = QLineEdit(self,"energyBox")
    self.energyBox.setReadOnly(1)

    layout3.addWidget(self.energyBox,0,0)

    self.textLabel2 = QLabel(self,"textLabel2")
    self.textLabel2.setAlignment(QLabel.AlignCenter)

    layout3.addWidget(self.textLabel2,1,1)

    self.textLabel1 = QLabel(self,"textLabel1")
    self.textLabel1.setAlignment(QLabel.AlignCenter)

    layout3.addWidget(self.textLabel1,1,0)
    AtomCorrespondenceDialogLayout.addLayout(layout3)

    layout1 = QHBoxLayout(None,0,6,"layout1")

    self.okButton = QPushButton(self,"okButton")
    layout1.addWidget(self.okButton)
    spacer1 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout1.addItem(spacer1)

    self.cancelButton = QPushButton(self,"cancelButton")
    layout1.addWidget(self.cancelButton)
    AtomCorrespondenceDialogLayout.addLayout(layout1)

    self.languageChange()

    self.resize(QSize(422,375).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.refineButton,SIGNAL("clicked()"),self.refineClicked)
    self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
    self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)
    self.connect(self.makeCorrespondenceButton,SIGNAL("clicked()"),self.makeCorrespondenceClicked)
    self.connect(self.showCorrespondenceBox,SIGNAL("stateChanged(int)"),self.showCorrespondenceChanged)
    self.connect(self.clearCorrespondenceButton,SIGNAL("clicked()"),self.clearCorrespondenceClicked)


  def languageChange(self):
    self.setCaption(self.__tr("Setup Atom-Atom Correspondences"))
    self.showCorrespondenceBox.setText(self.__tr("Show Correspondences?"))
    QToolTip.add(self.showCorrespondenceBox,self.__tr("Toggles display of the current correspondences in PyMol."))
    self.makeCorrespondenceButton.setText(self.__tr("Connect Picked Atoms"))
    QToolTip.add(self.makeCorrespondenceButton,self.__tr("Sets a correspondence between two picked atoms in PyMol."))
    self.clearCorrespondenceButton.setText(self.__tr("Clear Correspondences"))
    QToolTip.add(self.clearCorrespondenceButton,self.__tr("Clears all atom correspondences."))
    self.refineButton.setText(self.__tr("Refine Alignment"))
    QToolTip.add(self.refineButton,self.__tr("Refine the alignment using the existing atom-atom correspondence."))
    self.textLabel3.setText(self.__tr("RMSD"))
    self.rmsdBox.setInputMask(QString.null)
    QToolTip.add(self.rmsdBox,self.__tr("RMSD of the current alignment."))
    self.shapeScoreBox.setInputMask(QString.null)
    QToolTip.add(self.shapeScoreBox,self.__tr("Shape score for the current alignment."))
    self.energyBox.setInputMask(QString.null)
    QToolTip.add(self.energyBox,self.__tr("Relative energy of the current alignment."))
    self.textLabel2.setText(self.__tr("Shape Score"))
    self.textLabel1.setText(self.__tr("Relative Energy"))
    self.okButton.setText(self.__tr("OK"))
    QToolTip.add(self.okButton,self.__tr("Accept the current refined alignment and add it to the results pane."))
    self.cancelButton.setText(self.__tr("Cancel"))
    QToolTip.add(self.cancelButton,self.__tr("Cancel and discard the current refined alignment."))


  def refineClicked(self):
    print "AtomCorrespondenceDialog.refineClicked(): Not implemented yet"

  def showCorrespondenceChanged(self,a0):
    print "AtomCorrespondenceDialog.showCorrespondenceChanged(int): Not implemented yet"

  def makeCorrespondenceClicked(self):
    print "AtomCorrespondenceDialog.makeCorrespondenceClicked(): Not implemented yet"

  def clearCorrespondenceClicked(self):
    print "AtomCorrespondenceDialog.clearCorrespondenceClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("AtomCorrespondenceDialog",s,c)
