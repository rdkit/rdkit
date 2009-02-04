# Form implementation generated from reading ui file 'ReactionWindow.ui'
#
# Created: Fri Oct 25 14:41:57 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class ReactionWindow(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("reactionWindow")

    self.resize(663,508)
    self.setSizePolicy(QSizePolicy(5,5,0,0,self.sizePolicy().hasHeightForWidth()))
    self.setCaption(self.trUtf8("Reaction"))
    self.setSizeGripEnabled(0)

    reactionWindowLayout = QVBoxLayout(self,11,6,"reactionWindowLayout")

    Layout47 = QHBoxLayout(None,0,6,"Layout47")
    spacer = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout47.addItem(spacer)

    self.rxnPicture = QPushButton(self,"rxnPicture")
    self.rxnPicture.setSizePolicy(QSizePolicy(0,0,0,0,self.rxnPicture.sizePolicy().hasHeightForWidth()))
    self.rxnPicture.setMinimumSize(QSize(500,250))
    self.rxnPicture.setText(self.trUtf8("rxnPicture"))
    Layout47.addWidget(self.rxnPicture)
    spacer_2 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout47.addItem(spacer_2)
    reactionWindowLayout.addLayout(Layout47)

    self.infoFrame = QFrame(self,"infoFrame")
    self.infoFrame.setFrameShape(QFrame.StyledPanel)
    self.infoFrame.setFrameShadow(QFrame.Raised)
    infoFrameLayout = QVBoxLayout(self.infoFrame,11,6,"infoFrameLayout")

    Layout33 = QHBoxLayout(None,0,6,"Layout33")

    self.nameFormatBox = QLineEdit(self.infoFrame,"nameFormatBox")
    self.nameFormatBox.setText(self.trUtf8("%(react1)s_%(react2)s_%(productNum)d"))
    QToolTip.add(self.nameFormatBox,self.trUtf8("Python format string for product names."))
    Layout33.addWidget(self.nameFormatBox)

    self.prodLabel = QLabel(self.infoFrame,"prodLabel")
    self.prodLabel.setText(self.trUtf8("Product Name Format"))
    QToolTip.add(self.prodLabel,self.trUtf8("Python format string for product names."))
    Layout33.addWidget(self.prodLabel)
    infoFrameLayout.addLayout(Layout33)

    Layout34 = QHBoxLayout(None,0,6,"Layout34")

    self.multiProdCheck = QCheckBox(self.infoFrame,"multiProdCheck")
    self.multiProdCheck.setText(self.trUtf8("Multiple Products Allowed?"))
    QToolTip.add(self.multiProdCheck,self.trUtf8("Allow products to react again? This enables chain reactions."))
    Layout34.addWidget(self.multiProdCheck)
    spacer_3 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout34.addItem(spacer_3)
    infoFrameLayout.addLayout(Layout34)

    Layout35 = QHBoxLayout(None,0,6,"Layout35")

    self.maxCyclesSpin = QSpinBox(self.infoFrame,"maxCyclesSpin")
    self.maxCyclesSpin.setEnabled(1)
    self.maxCyclesSpin.setMinValue(1)
    QToolTip.add(self.maxCyclesSpin,self.trUtf8("maximum number of times the reaction will be allowed to cycle"))
    Layout35.addWidget(self.maxCyclesSpin)

    self.maxCyclesLabel = QLabel(self.infoFrame,"maxCyclesLabel")
    self.maxCyclesLabel.setEnabled(1)
    self.maxCyclesLabel.setText(self.trUtf8("Maximum cycles"))
    QToolTip.add(self.maxCyclesLabel,self.trUtf8("maximum number of times the reaction will be allowed to cycle"))
    Layout35.addWidget(self.maxCyclesLabel)

    self.intermediatesCheck = QCheckBox(self.infoFrame,"intermediatesCheck")
    self.intermediatesCheck.setEnabled(1)
    self.intermediatesCheck.setText(self.trUtf8("Retain Intermediates?"))
    QToolTip.add(self.intermediatesCheck,self.trUtf8("keep intermediates from chain reactions?"))
    Layout35.addWidget(self.intermediatesCheck)
    spacer_4 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout35.addItem(spacer_4)
    infoFrameLayout.addLayout(Layout35)

    self.Line1 = QFrame(self.infoFrame,"Line1")
    self.Line1.setProperty("frameShape",QVariant(QFrame.HLine))
    self.Line1.setFrameShadow(QFrame.Sunken)
    self.Line1.setFrameShape(QFrame.HLine)
    infoFrameLayout.addWidget(self.Line1)

    Layout38 = QHBoxLayout(None,0,6,"Layout38")

    self.nodeTextLabel = QLabel(self.infoFrame,"nodeTextLabel")
    self.nodeTextLabel.setText(self.trUtf8("Node Text"))
    QToolTip.add(self.nodeTextLabel,self.trUtf8("Text displayed on the canvas"))
    Layout38.addWidget(self.nodeTextLabel)

    self.nodeTextBox = QLineEdit(self.infoFrame,"nodeTextBox")
    self.nodeTextBox.setText(self.trUtf8("Rxn"))
    QToolTip.add(self.nodeTextBox,self.trUtf8("Text displayed on the canvas"))
    Layout38.addWidget(self.nodeTextBox)
    infoFrameLayout.addLayout(Layout38)
    reactionWindowLayout.addWidget(self.infoFrame)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer_5 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer_5)

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
    reactionWindowLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.rxnPicture,SIGNAL("clicked()"),self.rxnClick)
    self.connect(self.multiProdCheck,SIGNAL("pressed()"),self.multiStepClick)

  def rxnClick(self):
    print "ReactionWindow.rxnClick(): Not implemented yet"

  def multiStepClick(self):
    print "ReactionWindow.multiStepClick(): Not implemented yet"
