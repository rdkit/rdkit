# Form implementation generated from reading ui file 'FilterWindow.ui'
#
# Created: Fri Oct 25 15:09:45 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class FilterWindow(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("FilterWindow")

    self.resize(592,507)
    self.setCaption(self.trUtf8("Filter"))
    self.setSizeGripEnabled(0)
    QToolTip.add(self,self.trUtf8("num"))

    FilterWindowLayout = QVBoxLayout(self,11,6,"FilterWindowLayout")

    self.filterTab = QTabWidget(self,"filterTab")

    self.tab = QWidget(self.filterTab,"tab")

    LayoutWidget = QWidget(self.tab,"Layout41")
    LayoutWidget.setGeometry(QRect(11,321,544,54))
    Layout41 = QGridLayout(LayoutWidget,1,1,0,6,"Layout41")

    self.info_numInputsBox = QSpinBox(LayoutWidget,"info_numInputsBox")
    self.info_numInputsBox.setEnabled(0)
    QToolTip.add(self.info_numInputsBox,self.trUtf8("number of inputs to the node"))

    Layout41.addWidget(self.info_numInputsBox,0,0)

    self.info_inputLabel = QLabel(LayoutWidget,"info_inputLabel")
    self.info_inputLabel.setText(self.trUtf8("Number of Inputs"))
    QToolTip.add(self.info_inputLabel,self.trUtf8("number of inputs to the node"))

    Layout41.addWidget(self.info_inputLabel,0,1)

    self.info_numOutputsBox = QSpinBox(LayoutWidget,"info_numOutputsBox")
    self.info_numOutputsBox.setEnabled(0)
    QToolTip.add(self.info_numOutputsBox,self.trUtf8("number of outputs from the node"))

    Layout41.addWidget(self.info_numOutputsBox,1,0)

    self.info_outputLabel_2 = QLabel(LayoutWidget,"info_outputLabel_2")
    self.info_outputLabel_2.setText(self.trUtf8("Number of Outputs"))
    QToolTip.add(self.info_outputLabel_2,self.trUtf8("number of outputs from the node"))

    Layout41.addWidget(self.info_outputLabel_2,1,1)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout41.addItem(spacer,1,2)
    spacer_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout41.addItem(spacer_2,0,2)

    LayoutWidget_2 = QWidget(self.tab,"Layout66")
    LayoutWidget_2.setGeometry(QRect(11,11,544,116))
    Layout66 = QHBoxLayout(LayoutWidget_2,0,6,"Layout66")

    self.modeGroup = QButtonGroup(LayoutWidget_2,"modeGroup")
    self.modeGroup.setTitle(self.trUtf8("Mode"))

    self.codeButton = QRadioButton(self.modeGroup,"codeButton")
    self.codeButton.setGeometry(QRect(10,60,130,24))
    self.codeButton.setText(self.trUtf8("Code"))
    QToolTip.add(self.codeButton,self.trUtf8("Use a filter from a module"))

    self.moduleButton = QRadioButton(self.modeGroup,"moduleButton")
    self.moduleButton.setGeometry(QRect(10,30,130,24))
    self.moduleButton.setText(self.trUtf8("Module"))
    self.moduleButton.setChecked(1)
    QToolTip.add(self.moduleButton,self.trUtf8("Use a filter from a module"))
    Layout66.addWidget(self.modeGroup)
    spacer_3 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout66.addItem(spacer_3)

    self.moduleGroup = QGroupBox(self.tab,"moduleGroup")
    self.moduleGroup.setGeometry(QRect(11,133,544,150))
    self.moduleGroup.setTitle(self.trUtf8("Module"))

    LayoutWidget_3 = QWidget(self.moduleGroup,"Layout43")
    LayoutWidget_3.setGeometry(QRect(10,28,520,56))
    Layout43 = QGridLayout(LayoutWidget_3,1,1,0,6,"Layout43")

    Layout43_2 = QHBoxLayout(None,0,6,"Layout43_2")

    self.info_nameBox = QLineEdit(LayoutWidget_3,"info_nameBox")
    QToolTip.add(self.info_nameBox,self.trUtf8("file containing the filter function"))
    Layout43_2.addWidget(self.info_nameBox)

    self.info_nameButton = QToolButton(LayoutWidget_3,"info_nameButton")
    self.info_nameButton.setText(self.trUtf8("..."))
    Layout43_2.addWidget(self.info_nameButton)

    Layout43.addLayout(Layout43_2,0,1)

    self.info_Label_2_2_2 = QLabel(LayoutWidget_3,"info_Label_2_2_2")
    self.info_Label_2_2_2.setText(self.trUtf8("Method"))
    self.info_Label_2_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    QToolTip.add(self.info_Label_2_2_2,self.trUtf8("name of the method to call"))

    Layout43.addWidget(self.info_Label_2_2_2,1,0)

    self.info_methodBox = QLineEdit(LayoutWidget_3,"info_methodBox")
    self.info_methodBox.setEnabled(0)
    QToolTip.add(self.info_methodBox,self.trUtf8("name of the method to call"))

    Layout43.addWidget(self.info_methodBox,1,1)

    self.info_Label_2_3 = QLabel(LayoutWidget_3,"info_Label_2_3")
    self.info_Label_2_3.setText(self.trUtf8("File"))
    self.info_Label_2_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout43.addWidget(self.info_Label_2_3,0,0)
    self.filterTab.insertTab(self.tab,self.trUtf8("Info"))

    self.tab_2 = QWidget(self.filterTab,"tab_2")
    tabLayout = QVBoxLayout(self.tab_2,11,6,"tabLayout")

    self.code_editBox = QTextEdit(self.tab_2,"code_editBox")
    self.code_editBox.setText(self.trUtf8("def filter(self,inputValues):\n"
"  \"\"\" \n"
"   inputValues is a python list with the values of this node's inputs  \n"
"   the function should return an integer\n"
"  \"\"\"\n"
"  sum = 0\n"
"  for val in inputValues:\n"
"     sum += int(val)\n"
"  if sum > 2:\n"
"     return 1\n"
"  else:\n"
"     return 0\n"
""))
    QToolTip.add(self.code_editBox,self.trUtf8("enter python code for the filter here"))
    tabLayout.addWidget(self.code_editBox)
    self.filterTab.insertTab(self.tab_2,self.trUtf8("Code"))
    FilterWindowLayout.addWidget(self.filterTab)

    Layout40 = QHBoxLayout(None,0,6,"Layout40")

    self.nodeTextLabel = QLabel(self,"nodeTextLabel")
    self.nodeTextLabel.setText(self.trUtf8("Node Text"))
    QToolTip.add(self.nodeTextLabel,self.trUtf8("Text displayed on the canvas"))
    Layout40.addWidget(self.nodeTextLabel)

    self.nodeTextBox = QLineEdit(self,"nodeTextBox")
    self.nodeTextBox.setText(self.trUtf8("Filter"))
    QToolTip.add(self.nodeTextBox,self.trUtf8("Text displayed on the canvas"))
    Layout40.addWidget(self.nodeTextBox)
    FilterWindowLayout.addLayout(Layout40)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer_4 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer_4)

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
    FilterWindowLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.info_nameButton,SIGNAL("clicked()"),self.fileClick)

  def fileClick(self):
    print "FilterWindow.fileClick(): Not implemented yet"

  def modeModuleClick(self):
    print "FilterWindow.modeModuleClick(): Not implemented yet"

  def modeCodeClick(self):
    print "FilterWindow.modeCodeClick(): Not implemented yet"
