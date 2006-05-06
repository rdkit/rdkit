# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MolEditDisplayWidget2.ui'
#
# Created: Tue Jan 24 14:43:33 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.11
#
# WARNING! All changes made in this file will be lost!


from qt import *


class MolEditDisplayWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if not name:
      self.setName("MolEditDisplayWidget")

    self.setAcceptDrops(1)

    MolEditDisplayWidgetLayout = QVBoxLayout(self,2,2,"MolEditDisplayWidgetLayout")

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.fromCDXButton = QPushButton(self,"fromCDXButton")
    self.fromCDXButton.setEnabled(0)
    Layout1.addWidget(self.fromCDXButton)
    Spacer1 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(Spacer1)

    self.toCDXButton = QPushButton(self,"toCDXButton")
    self.toCDXButton.setEnabled(0)
    Layout1.addWidget(self.toCDXButton)
    MolEditDisplayWidgetLayout.addLayout(Layout1)

    self.formatGroup = QButtonGroup(self,"formatGroup")
    self.formatGroup.setSizePolicy(QSizePolicy(5,0,0,0,self.formatGroup.sizePolicy().hasHeightForWidth()))
    self.formatGroup.setColumnLayout(0,Qt.Vertical)
    self.formatGroup.layout().setSpacing(6)
    self.formatGroup.layout().setMargin(6)
    formatGroupLayout = QVBoxLayout(self.formatGroup.layout())
    formatGroupLayout.setAlignment(Qt.AlignTop)

    Layout33 = QGridLayout(None,1,1,0,6,"Layout33")

    self.placeholder1 = QLabel(self.formatGroup,"placeholder1")

    Layout33.addWidget(self.placeholder1,1,2)

    self.fileButton = QToolButton(self.formatGroup,"fileButton")

    Layout33.addWidget(self.fileButton,0,2)

    self.fileRadio = QRadioButton(self.formatGroup,"fileRadio")
    self.fileRadio.setChecked(1)

    Layout33.addWidget(self.fileRadio,0,0)

    self.fileEdit = QLineEdit(self.formatGroup,"fileEdit")
    self.fileEdit.setAcceptDrops(0)
    self.fileEdit.setDragEnabled(1)

    Layout33.addWidget(self.fileEdit,0,1)

    self.smilesRadio = QRadioButton(self.formatGroup,"smilesRadio")

    Layout33.addWidget(self.smilesRadio,1,0)

    self.smilesEdit = QLineEdit(self.formatGroup,"smilesEdit")
    self.smilesEdit.setAcceptDrops(0)
    self.smilesEdit.setDragEnabled(1)

    Layout33.addWidget(self.smilesEdit,1,1)
    formatGroupLayout.addLayout(Layout33)
    MolEditDisplayWidgetLayout.addWidget(self.formatGroup)

    self.languageChange()

    self.resize(QSize(454,411).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.fromCDXButton,SIGNAL("clicked()"),self.fromCDXClick)
    self.connect(self.toCDXButton,SIGNAL("clicked()"),self.toCDXClick)
    self.connect(self.smilesEdit,SIGNAL("textChanged(const QString&)"),self.smilesEdited)
    self.connect(self.fileButton,SIGNAL("clicked()"),self.fileButtonClicked)


  def languageChange(self):
    self.setCaption(self.__tr("MolEditDisplay"))
    self.fromCDXButton.setText(self.__tr("From ChemDraw"))
    QToolTip.add(self.fromCDXButton,self.__tr("Grab the contents of the active ChemDraw document."))
    self.toCDXButton.setText(self.__tr("To ChemDraw"))
    QToolTip.add(self.toCDXButton,self.__tr("Open the current molecule in ChemDraw for editing."))
    self.formatGroup.setTitle(self.__tr("Input"))
    self.placeholder1.setText(QString.null)
    self.fileButton.setText(self.__tr("..."))
    QToolTip.add(self.fileButton,self.__tr("Browse to the input file."))
    self.fileRadio.setText(self.__tr("File"))
    QToolTip.add(self.fileEdit,self.__tr("Enter the name of the file containing a molecule.  Can be either a .mol or .mol2 file."))
    self.smilesRadio.setText(self.__tr("SMILES"))
    QToolTip.add(self.smilesEdit,self.__tr("Enter the Daylight SMILES for the molecule."))


  def fromCDXClick(self):
    print "MolEditDisplayWidget.fromCDXClick(): Not implemented yet"

  def toCDXClick(self):
    print "MolEditDisplayWidget.toCDXClick(): Not implemented yet"

  def smilesEdited(self):
    print "MolEditDisplayWidget.smilesEdited(): Not implemented yet"

  def fileButtonClicked(self):
    print "MolEditDisplayWidget.fileButtonClicked(): Not implemented yet"

  def nextButtonClicked(self):
    print "MolEditDisplayWidget.nextButtonClicked(): Not implemented yet"

  def previousButtonClicked(self):
    print "MolEditDisplayWidget.previousButtonClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("MolEditDisplayWidget",s,c)
