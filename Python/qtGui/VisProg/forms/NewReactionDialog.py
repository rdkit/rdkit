# Form implementation generated from reading ui file 'NewReactionDialog.ui'
#
# Created: Fri Oct 25 15:52:10 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class NewReactionDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("NewReactionDialog")

    self.resize(516,259)
    self.setCaption(self.trUtf8("New Reaction"))
    self.setSizeGripEnabled(0)


    LayoutWidget = QWidget(self,"Layout1")
    LayoutWidget.setGeometry(QRect(10,220,500,38))
    Layout1 = QHBoxLayout(LayoutWidget,0,6,"Layout1")

    self.buttonHelp = QPushButton(LayoutWidget,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer)

    self.buttonOk = QPushButton(LayoutWidget,"buttonOk")
    self.buttonOk.setText(self.trUtf8("OK"))
    self.buttonOk.setAccel(0)
    self.buttonOk.setAutoDefault(1)
    self.buttonOk.setDefault(1)
    Layout1.addWidget(self.buttonOk)

    self.buttonCancel = QPushButton(LayoutWidget,"buttonCancel")
    self.buttonCancel.setText(self.trUtf8("Cancel"))
    self.buttonCancel.setAccel(0)
    self.buttonCancel.setAutoDefault(1)
    Layout1.addWidget(self.buttonCancel)

    LayoutWidget_2 = QWidget(self,"Layout57")
    LayoutWidget_2.setGeometry(QRect(31,11,480,200))
    Layout57 = QVBoxLayout(LayoutWidget_2,0,6,"Layout57")

    Layout50 = QGridLayout(None,1,1,0,6,"Layout50")

    self.imgFileBrowse = QToolButton(LayoutWidget_2,"imgFileBrowse")
    self.imgFileBrowse.setText(self.trUtf8("..."))
    QToolTip.add(self.imgFileBrowse,self.trUtf8("open file dialog"))

    Layout50.addWidget(self.imgFileBrowse,1,2)

    self.imgLabel = QLabel(LayoutWidget_2,"imgLabel")
    self.imgLabel.setText(self.trUtf8("Image File"))
    self.imgLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout50.addWidget(self.imgLabel,1,0)

    self.imgField = QLineEdit(LayoutWidget_2,"imgField")
    QToolTip.add(self.imgField,self.trUtf8("name of the image file to be used"))

    Layout50.addWidget(self.imgField,1,1)

    self.rxnFileBrowse = QToolButton(LayoutWidget_2,"rxnFileBrowse")
    self.rxnFileBrowse.setText(self.trUtf8("..."))
    QToolTip.add(self.rxnFileBrowse,self.trUtf8("open file dialog"))

    Layout50.addWidget(self.rxnFileBrowse,0,2)

    self.rxnLabel = QLabel(LayoutWidget_2,"rxnLabel")
    self.rxnLabel.setText(self.trUtf8("Reaction File"))
    self.rxnLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout50.addWidget(self.rxnLabel,0,0)

    self.rxnField = QLineEdit(LayoutWidget_2,"rxnField")
    QToolTip.add(self.rxnField,self.trUtf8("name of the reaction file to be used"))

    Layout50.addWidget(self.rxnField,0,1)
    Layout57.addLayout(Layout50)

    self.Line2_2_2 = QFrame(LayoutWidget_2,"Line2_2_2")
    self.Line2_2_2.setFrameShape(QFrame.HLine)
    self.Line2_2_2.setFrameShadow(QFrame.Sunken)
    self.Line2_2_2.setFrameShape(QFrame.HLine)
    Layout57.addWidget(self.Line2_2_2)

    Layout51 = QHBoxLayout(None,0,6,"Layout51")

    self.smirksLabel = QLabel(LayoutWidget_2,"smirksLabel")
    self.smirksLabel.setText(self.trUtf8("SMIRKS"))
    self.smirksLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    Layout51.addWidget(self.smirksLabel)

    self.smirksField = QLineEdit(LayoutWidget_2,"smirksField")
    QToolTip.add(self.smirksField,self.trUtf8("enter SMIRKS for the reaction here"))
    Layout51.addWidget(self.smirksField)
    Layout57.addLayout(Layout51)

    self.Line2_2 = QFrame(LayoutWidget_2,"Line2_2")
    self.Line2_2.setFrameShape(QFrame.HLine)
    self.Line2_2.setFrameShadow(QFrame.Sunken)
    self.Line2_2.setFrameShape(QFrame.HLine)
    Layout57.addWidget(self.Line2_2)

    Layout56 = QGridLayout(None,1,1,0,6,"Layout56")

    self.prodLabel = QLabel(LayoutWidget_2,"prodLabel")
    self.prodLabel.setText(self.trUtf8("Num Products"))
    self.prodLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout56.addWidget(self.prodLabel,1,0)

    self.reactLabel = QLabel(LayoutWidget_2,"reactLabel")
    self.reactLabel.setText(self.trUtf8("Num Reactants"))
    self.reactLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout56.addWidget(self.reactLabel,0,0)
    spacer_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout56.addItem(spacer_2,1,2)

    self.nProdSpinner = QSpinBox(LayoutWidget_2,"nProdSpinner")
    QToolTip.add(self.nProdSpinner,self.trUtf8("number of reaction products"))

    Layout56.addWidget(self.nProdSpinner,1,1)

    self.nReactSpinner = QSpinBox(LayoutWidget_2,"nReactSpinner")
    QToolTip.add(self.nReactSpinner,self.trUtf8("number of reactants"))

    Layout56.addWidget(self.nReactSpinner,0,1)
    spacer_3 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout56.addItem(spacer_3,0,2)
    Layout57.addLayout(Layout56)

    self.radio = QButtonGroup(self,"radio")
    self.radio.setGeometry(QRect(0,0,30,220))
    self.radio.setLineWidth(0)
    self.radio.setTitle(self.trUtf8(""))
    self.radio.setExclusive(1)

    self.fileRadio = QRadioButton(self.radio,"fileRadio")
    self.fileRadio.setGeometry(QRect(10,41,20,16))
    self.fileRadio.setText(self.trUtf8(""))
    self.fileRadio.setChecked(1)

    self.smirksRadio = QRadioButton(self.radio,"smirksRadio")
    self.smirksRadio.setGeometry(QRect(10,109,20,16))
    self.smirksRadio.setText(self.trUtf8(""))

    self.numberRadio = QRadioButton(self.radio,"numberRadio")
    self.numberRadio.setGeometry(QRect(10,176,20,16))
    self.numberRadio.setText(self.trUtf8(""))

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.imgFileBrowse,SIGNAL("clicked()"),self.imgFileClick)
    self.connect(self.nReactSpinner,SIGNAL("valueChanged(int)"),self.spinClick)
    self.connect(self.nProdSpinner,SIGNAL("valueChanged(int)"),self.spinClick)
    self.connect(self.rxnFileBrowse,SIGNAL("clicked()"),self.rxnFileClick)

    self.setTabOrder(self.rxnField,self.imgField)
    self.setTabOrder(self.imgField,self.smirksField)
    self.setTabOrder(self.smirksField,self.nReactSpinner)
    self.setTabOrder(self.nReactSpinner,self.nProdSpinner)
    self.setTabOrder(self.nProdSpinner,self.buttonOk)
    self.setTabOrder(self.buttonOk,self.buttonCancel)
    self.setTabOrder(self.buttonCancel,self.buttonHelp)
    self.setTabOrder(self.buttonHelp,self.numberRadio)
    self.setTabOrder(self.numberRadio,self.smirksRadio)
    self.setTabOrder(self.smirksRadio,self.fileRadio)

  def imgFileClick(self):
    print "NewReactionDialog.imgFileClick(): Not implemented yet"

  def rxnFileClick(self):
    print "NewReactionDialog.rxnFileClick(): Not implemented yet"

  def spinClick(self):
    print "NewReactionDialog.spinClick(): Not implemented yet"
