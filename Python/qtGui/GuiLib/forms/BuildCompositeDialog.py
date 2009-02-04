# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/BuildCompositeDialog.ui'
#
# Created: Sun Apr 23 08:10:30 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.11
#
# WARNING! All changes made in this file will be lost!


from qt import *
from qttable import QTable

image0_data = [
"16 16 6 1",
". c None",
"# c #1414ff",
"a c #143cff",
"b c #6464ff",
"c c #8c8cff",
"d c #b4b4ff",
"................",
".....######.....",
"....#aaaaaa#....",
"...#abbbbbba#...",
"..#abbbbbbbba#..",
".#abbbccccbbba#.",
".#abbccccccbba#.",
".#abbccddccbba#.",
".#abbccddccbba#.",
".#abbccccccbba#.",
".#abbbccccbbba#.",
"..#abbbbbbbba#..",
"...#abbbbbba#...",
"....#aaaaaa#....",
".....######.....",
"................"
]

class CompositeBuilderDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("CompositeBuilder")

    self.setSizePolicy(QSizePolicy(5,3,0,0,self.sizePolicy().hasHeightForWidth()))
    self.setIcon(self.image0)
    self.setSizeGripEnabled(0)

    CompositeBuilderLayout = QVBoxLayout(self,4,2,"CompositeBuilderLayout")

    self.mainTab = QTabWidget(self,"mainTab")
    self.mainTab.setEnabled(1)

    self.paramPage = QWidget(self.mainTab,"paramPage")
    paramPageLayout = QVBoxLayout(self.paramPage,4,2,"paramPageLayout")

    layout33 = QHBoxLayout(None,0,2,"layout33")

    layout94 = QVBoxLayout(None,0,2,"layout94")

    Layout23 = QHBoxLayout(None,0,6,"Layout23")

    self.params_count = QSpinBox(self.paramPage,"params_count")
    self.params_count.setSizePolicy(QSizePolicy(0,0,0,0,self.params_count.sizePolicy().hasHeightForWidth()))
    self.params_count.setMaxValue(1000)
    self.params_count.setLineStep(5)
    self.params_count.setValue(10)
    Layout23.addWidget(self.params_count)

    self.TextLabel4 = QLabel(self.paramPage,"TextLabel4")
    Layout23.addWidget(self.TextLabel4)
    layout94.addLayout(Layout23)

    self.params_bayes = QCheckBox(self.paramPage,"params_bayes")
    self.params_bayes.setChecked(0)
    layout94.addWidget(self.params_bayes)

    self.GroupBox1 = QGroupBox(self.paramPage,"GroupBox1")
    self.GroupBox1.setColumnLayout(0,Qt.Vertical)
    self.GroupBox1.layout().setSpacing(4)
    self.GroupBox1.layout().setMargin(2)
    GroupBox1Layout = QVBoxLayout(self.GroupBox1.layout())
    GroupBox1Layout.setAlignment(Qt.AlignTop)

    Layout19 = QHBoxLayout(None,0,6,"Layout19")

    self.params_split = QCheckBox(self.GroupBox1,"params_split")
    self.params_split.setChecked(1)
    Layout19.addWidget(self.params_split)

    self.params_splitFrac = QLineEdit(self.GroupBox1,"params_splitFrac")
    self.params_splitFrac.setSizePolicy(QSizePolicy(0,0,0,0,self.params_splitFrac.sizePolicy().hasHeightForWidth()))
    self.params_splitFrac.setMaximumSize(QSize(40,32767))
    Layout19.addWidget(self.params_splitFrac)

    self.TextLabel4_2_2_2 = QLabel(self.GroupBox1,"TextLabel4_2_2_2")
    Layout19.addWidget(self.TextLabel4_2_2_2)
    GroupBox1Layout.addLayout(Layout19)
    layout94.addWidget(self.GroupBox1)

    self.ButtonGroup1 = QButtonGroup(self.paramPage,"ButtonGroup1")
    self.ButtonGroup1.setColumnLayout(0,Qt.Vertical)
    self.ButtonGroup1.layout().setSpacing(4)
    self.ButtonGroup1.layout().setMargin(2)
    ButtonGroup1Layout = QVBoxLayout(self.ButtonGroup1.layout())
    ButtonGroup1Layout.setAlignment(Qt.AlignTop)

    self.params_dataNone = QRadioButton(self.ButtonGroup1,"params_dataNone")
    self.params_dataNone.setChecked(1)
    ButtonGroup1Layout.addWidget(self.params_dataNone)

    self.params_dataShuffle = QRadioButton(self.ButtonGroup1,"params_dataShuffle")
    ButtonGroup1Layout.addWidget(self.params_dataShuffle)

    self.params_dataRandomize = QRadioButton(self.ButtonGroup1,"params_dataRandomize")
    ButtonGroup1Layout.addWidget(self.params_dataRandomize)
    layout94.addWidget(self.ButtonGroup1)
    Spacer16_2 = QSpacerItem(20,182,QSizePolicy.Minimum,QSizePolicy.Expanding)
    layout94.addItem(Spacer16_2)

    Layout171 = QHBoxLayout(None,0,6,"Layout171")

    self.params_actBounds = QLineEdit(self.paramPage,"params_actBounds")
    self.params_actBounds.setSizePolicy(QSizePolicy(0,0,0,0,self.params_actBounds.sizePolicy().hasHeightForWidth()))
    self.params_actBounds.setMaximumSize(QSize(40,32767))
    Layout171.addWidget(self.params_actBounds)

    self.TextLabel4_2_2_2_2 = QLabel(self.paramPage,"TextLabel4_2_2_2_2")
    self.TextLabel4_2_2_2_2.setMaximumSize(QSize(32767,32767))
    Layout171.addWidget(self.TextLabel4_2_2_2_2)
    layout94.addLayout(Layout171)

    Layout17 = QHBoxLayout(None,0,6,"Layout17")

    self.params_threshold = QLineEdit(self.paramPage,"params_threshold")
    self.params_threshold.setEnabled(1)
    self.params_threshold.setSizePolicy(QSizePolicy(0,0,0,0,self.params_threshold.sizePolicy().hasHeightForWidth()))
    self.params_threshold.setMaximumSize(QSize(40,32767))
    Layout17.addWidget(self.params_threshold)

    self.TextLabel4_2_2 = QLabel(self.paramPage,"TextLabel4_2_2")
    self.TextLabel4_2_2.setMaximumSize(QSize(32767,32767))
    Layout17.addWidget(self.TextLabel4_2_2)
    layout94.addLayout(Layout17)
    layout33.addLayout(layout94)
    Spacer14_2 = QSpacerItem(159,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout33.addItem(Spacer14_2)

    layout32 = QVBoxLayout(None,0,2,"layout32")

    self.GroupBox4 = QGroupBox(self.paramPage,"GroupBox4")
    self.GroupBox4.setColumnLayout(0,Qt.Vertical)
    self.GroupBox4.layout().setSpacing(6)
    self.GroupBox4.layout().setMargin(11)
    GroupBox4Layout = QVBoxLayout(self.GroupBox4.layout())
    GroupBox4Layout.setAlignment(Qt.AlignTop)

    self.params_multiplebuilds = QCheckBox(self.GroupBox4,"params_multiplebuilds")
    GroupBox4Layout.addWidget(self.params_multiplebuilds)

    Layout12 = QHBoxLayout(None,0,6,"Layout12")
    Spacer12 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout12.addItem(Spacer12)

    self.params_multiplebuildcount = QSpinBox(self.GroupBox4,"params_multiplebuildcount")
    self.params_multiplebuildcount.setEnabled(0)
    self.params_multiplebuildcount.setSizePolicy(QSizePolicy(0,0,0,0,self.params_multiplebuildcount.sizePolicy().hasHeightForWidth()))
    self.params_multiplebuildcount.setMaxValue(1000)
    self.params_multiplebuildcount.setMinValue(1)
    self.params_multiplebuildcount.setValue(10)
    Layout12.addWidget(self.params_multiplebuildcount)

    self.TextLabel4_2_3 = QLabel(self.GroupBox4,"TextLabel4_2_3")
    Layout12.addWidget(self.TextLabel4_2_3)
    GroupBox4Layout.addLayout(Layout12)
    layout32.addWidget(self.GroupBox4)

    self.GroupBox2 = QGroupBox(self.paramPage,"GroupBox2")
    self.GroupBox2.setColumnLayout(0,Qt.Vertical)
    self.GroupBox2.layout().setSpacing(4)
    self.GroupBox2.layout().setMargin(2)
    GroupBox2Layout = QVBoxLayout(self.GroupBox2.layout())
    GroupBox2Layout.setAlignment(Qt.AlignTop)

    self.params_filter = QCheckBox(self.GroupBox2,"params_filter")
    GroupBox2Layout.addWidget(self.params_filter)

    layout31 = QGridLayout(None,1,1,0,2,"layout31")

    self.params_filterVal = QSpinBox(self.GroupBox2,"params_filterVal")
    self.params_filterVal.setEnabled(0)

    layout31.addWidget(self.params_filterVal,0,1)
    Spacer11 = QSpacerItem(52,16,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout31.addItem(Spacer11,1,0)

    self.params_filterModels = QCheckBox(self.GroupBox2,"params_filterModels")
    self.params_filterModels.setEnabled(0)
    self.params_filterModels.setChecked(1)

    layout31.addMultiCellWidget(self.params_filterModels,2,2,1,2)
    Spacer10 = QSpacerItem(52,16,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout31.addItem(Spacer10,0,0)

    self.TextLabel4_2_2_2_2_2 = QLabel(self.GroupBox2,"TextLabel4_2_2_2_2_2")

    layout31.addMultiCellWidget(self.TextLabel4_2_2_2_2_2,1,1,2,3)
    spacer21_2 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout31.addItem(spacer21_2,2,0)

    self.TextLabel1 = QLabel(self.GroupBox2,"TextLabel1")

    layout31.addMultiCellWidget(self.TextLabel1,0,0,2,3)

    self.params_filterFrac = QLineEdit(self.GroupBox2,"params_filterFrac")
    self.params_filterFrac.setEnabled(0)
    self.params_filterFrac.setSizePolicy(QSizePolicy(1,0,0,0,self.params_filterFrac.sizePolicy().hasHeightForWidth()))
    self.params_filterFrac.setMaximumSize(QSize(40,32767))

    layout31.addWidget(self.params_filterFrac,1,1)
    GroupBox2Layout.addLayout(layout31)
    layout32.addWidget(self.GroupBox2)

    self.GroupBox3 = QGroupBox(self.paramPage,"GroupBox3")
    self.GroupBox3.setColumnLayout(0,Qt.Vertical)
    self.GroupBox3.layout().setSpacing(4)
    self.GroupBox3.layout().setMargin(2)
    GroupBox3Layout = QVBoxLayout(self.GroupBox3.layout())
    GroupBox3Layout.setAlignment(Qt.AlignTop)

    Layout21 = QVBoxLayout(None,0,6,"Layout21")

    self.params_lock = QCheckBox(self.GroupBox3,"params_lock")
    Layout21.addWidget(self.params_lock)

    Layout20 = QGridLayout(None,1,1,0,6,"Layout20")

    self.params_lockV2 = QLineEdit(self.GroupBox3,"params_lockV2")
    self.params_lockV2.setEnabled(0)
    self.params_lockV2.setMaximumSize(QSize(40,32767))

    Layout20.addWidget(self.params_lockV2,1,1)
    Spacer9 = QSpacerItem(56,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer9,1,0)

    self.TextLabel1_2 = QLabel(self.GroupBox3,"TextLabel1_2")

    Layout20.addWidget(self.TextLabel1_2,0,2)

    self.TextLabel1_2_2 = QLabel(self.GroupBox3,"TextLabel1_2_2")

    Layout20.addWidget(self.TextLabel1_2_2,1,2)

    self.params_lockV1 = QLineEdit(self.GroupBox3,"params_lockV1")
    self.params_lockV1.setEnabled(0)
    self.params_lockV1.setMaximumSize(QSize(40,32767))

    Layout20.addWidget(self.params_lockV1,0,1)
    Spacer8 = QSpacerItem(56,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer8,0,0)
    Layout21.addLayout(Layout20)
    GroupBox3Layout.addLayout(Layout21)
    layout32.addWidget(self.GroupBox3)
    Spacer17 = QSpacerItem(20,20,QSizePolicy.Minimum,QSizePolicy.Expanding)
    layout32.addItem(Spacer17)
    layout33.addLayout(layout32)
    paramPageLayout.addLayout(layout33)

    Layout31 = QGridLayout(None,1,1,0,6,"Layout31")

    self.TextLabel1_3 = QLabel(self.paramPage,"TextLabel1_3")
    self.TextLabel1_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout31.addWidget(self.TextLabel1_3,0,0)

    self.TextLabel1_3_2 = QLabel(self.paramPage,"TextLabel1_3_2")
    self.TextLabel1_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout31.addWidget(self.TextLabel1_3_2,1,0)

    self.params_note = QLineEdit(self.paramPage,"params_note")

    Layout31.addWidget(self.params_note,1,1)

    self.params_persisttable = QLineEdit(self.paramPage,"params_persisttable")

    Layout31.addWidget(self.params_persisttable,0,1)
    paramPageLayout.addLayout(Layout31)
    self.mainTab.insertTab(self.paramPage,QString(""))

    self.TabPage = QWidget(self.mainTab,"TabPage")
    TabPageLayout = QHBoxLayout(self.TabPage,4,2,"TabPageLayout")

    self.buttonGroup2 = QButtonGroup(self.TabPage,"buttonGroup2")
    self.buttonGroup2.setColumnLayout(0,Qt.Vertical)
    self.buttonGroup2.layout().setSpacing(2)
    self.buttonGroup2.layout().setMargin(4)
    buttonGroup2Layout = QHBoxLayout(self.buttonGroup2.layout())
    buttonGroup2Layout.setAlignment(Qt.AlignTop)

    layout43 = QVBoxLayout(None,0,2,"layout43")

    self.params_treeRadio = QRadioButton(self.buttonGroup2,"params_treeRadio")
    self.params_treeRadio.setChecked(1)
    layout43.addWidget(self.params_treeRadio)

    layout42 = QHBoxLayout(None,0,2,"layout42")
    spacer16 = QSpacerItem(79,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout42.addItem(spacer16)

    layout26 = QVBoxLayout(None,0,2,"layout26")

    Layout22 = QHBoxLayout(None,0,6,"Layout22")

    self.params_depth = QSpinBox(self.buttonGroup2,"params_depth")
    self.params_depth.setSizePolicy(QSizePolicy(0,0,0,0,self.params_depth.sizePolicy().hasHeightForWidth()))
    self.params_depth.setMaxValue(1000)
    self.params_depth.setMinValue(-1)
    self.params_depth.setValue(-1)
    Layout22.addWidget(self.params_depth)

    self.TextLabel4_2 = QLabel(self.buttonGroup2,"TextLabel4_2")
    Layout22.addWidget(self.TextLabel4_2)
    layout26.addLayout(Layout22)

    self.params_autoBounds = QCheckBox(self.buttonGroup2,"params_autoBounds")
    self.params_autoBounds.setChecked(1)
    layout26.addWidget(self.params_autoBounds)

    self.params_greedy = QCheckBox(self.buttonGroup2,"params_greedy")
    layout26.addWidget(self.params_greedy)

    self.params_recycle = QCheckBox(self.buttonGroup2,"params_recycle")
    self.params_recycle.setChecked(0)
    layout26.addWidget(self.params_recycle)

    self.params_prune = QCheckBox(self.buttonGroup2,"params_prune")
    self.params_prune.setChecked(0)
    layout26.addWidget(self.params_prune)

    layout25 = QHBoxLayout(None,0,2,"layout25")

    self.params_randomDescriptors = QSpinBox(self.buttonGroup2,"params_randomDescriptors")
    self.params_randomDescriptors.setSizePolicy(QSizePolicy(0,0,0,0,self.params_randomDescriptors.sizePolicy().hasHeightForWidth()))
    self.params_randomDescriptors.setMaxValue(1000)
    self.params_randomDescriptors.setMinValue(0)
    self.params_randomDescriptors.setValue(0)
    layout25.addWidget(self.params_randomDescriptors)

    self.TextLabel4_2_4 = QLabel(self.buttonGroup2,"TextLabel4_2_4")
    layout25.addWidget(self.TextLabel4_2_4)
    layout26.addLayout(layout25)
    layout42.addLayout(layout26)
    layout43.addLayout(layout42)

    self.params_knnRadio = QRadioButton(self.buttonGroup2,"params_knnRadio")
    self.params_knnRadio.setEnabled(1)
    layout43.addWidget(self.params_knnRadio)

    layout39 = QHBoxLayout(None,0,2,"layout39")
    spacer18 = QSpacerItem(30,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout39.addItem(spacer18)

    layout24 = QVBoxLayout(None,0,2,"layout24")

    layout19 = QHBoxLayout(None,0,2,"layout19")

    self.params_knnKVal = QSpinBox(self.buttonGroup2,"params_knnKVal")
    self.params_knnKVal.setEnabled(0)
    self.params_knnKVal.setMinValue(1)
    self.params_knnKVal.setValue(5)
    layout19.addWidget(self.params_knnKVal)

    self.textLabel1 = QLabel(self.buttonGroup2,"textLabel1")
    self.textLabel1.setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
    layout19.addWidget(self.textLabel1)
    layout24.addLayout(layout19)

    layout20 = QHBoxLayout(None,0,2,"layout20")

    self.params_knnMetric = QComboBox(0,self.buttonGroup2,"params_knnMetric")
    self.params_knnMetric.setEnabled(0)
    layout20.addWidget(self.params_knnMetric)

    self.textLabel1_2 = QLabel(self.buttonGroup2,"textLabel1_2")
    self.textLabel1_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
    layout20.addWidget(self.textLabel1_2)
    layout24.addLayout(layout20)
    layout39.addLayout(layout24)
    layout43.addLayout(layout39)

    self.params_naiveBayesRadio = QRadioButton(self.buttonGroup2,"params_naiveBayesRadio")
    self.params_naiveBayesRadio.setEnabled(1)
    layout43.addWidget(self.params_naiveBayesRadio)

    layout40 = QHBoxLayout(None,0,2,"layout40")
    spacer18_2 = QSpacerItem(30,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout40.addItem(spacer18_2)

    layout34 = QHBoxLayout(None,0,2,"layout34")

    self.params_naiveBayesM = QLineEdit(self.buttonGroup2,"params_naiveBayesM")
    self.params_naiveBayesM.setEnabled(0)
    layout34.addWidget(self.params_naiveBayesM)

    self.textLabel4_3 = QLabel(self.buttonGroup2,"textLabel4_3")
    layout34.addWidget(self.textLabel4_3)
    layout40.addLayout(layout34)
    layout43.addLayout(layout40)
    spacer21 = QSpacerItem(20,43,QSizePolicy.Minimum,QSizePolicy.Expanding)
    layout43.addItem(spacer21)
    buttonGroup2Layout.addLayout(layout43)
    spacer59 = QSpacerItem(50,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    buttonGroup2Layout.addItem(spacer59)

    layout93 = QVBoxLayout(None,0,2,"layout93")

    layout92 = QVBoxLayout(None,0,2,"layout92")

    self.params_svmRadio = QRadioButton(self.buttonGroup2,"params_svmRadio")
    self.params_svmRadio.setEnabled(0)
    layout92.addWidget(self.params_svmRadio)

    layout23 = QHBoxLayout(None,0,2,"layout23")
    spacer17 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout23.addItem(spacer17)

    layout22 = QGridLayout(None,1,1,0,2,"layout22")

    self.params_svmEps = QLineEdit(self.buttonGroup2,"params_svmEps")
    self.params_svmEps.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmEps,4,4,0,1)

    self.params_svmCost = QLineEdit(self.buttonGroup2,"params_svmCost")
    self.params_svmCost.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmCost,8,8,0,1)

    self.textLabel1_3 = QLabel(self.buttonGroup2,"textLabel1_3")
    self.textLabel1_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)

    layout22.addMultiCellWidget(self.textLabel1_3,3,3,1,2)

    self.textLabel4_2_2 = QLabel(self.buttonGroup2,"textLabel4_2_2")
    self.textLabel4_2_2.setAlignment(QLabel.AlignVCenter)

    layout22.addWidget(self.textLabel4_2_2,6,2)

    self.params_svmGamma = QLineEdit(self.buttonGroup2,"params_svmGamma")
    self.params_svmGamma.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmGamma,7,7,0,1)

    self.params_svmCoeff = QLineEdit(self.buttonGroup2,"params_svmCoeff")
    self.params_svmCoeff.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmCoeff,6,6,0,1)

    self.textLabel6 = QLabel(self.buttonGroup2,"textLabel6")

    layout22.addWidget(self.textLabel6,8,2)

    self.textLabel5 = QLabel(self.buttonGroup2,"textLabel5")
    self.textLabel5.setAlignment(QLabel.AlignVCenter)

    layout22.addWidget(self.textLabel5,7,2)

    self.params_svmType = QComboBox(0,self.buttonGroup2,"params_svmType")
    self.params_svmType.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmType,1,1,0,1)

    self.textLabel4 = QLabel(self.buttonGroup2,"textLabel4")

    layout22.addWidget(self.textLabel4,4,2)

    self.textLabel2_2 = QLabel(self.buttonGroup2,"textLabel2_2")

    layout22.addWidget(self.textLabel2_2,1,2)

    self.textLabel2 = QLabel(self.buttonGroup2,"textLabel2")

    layout22.addWidget(self.textLabel2,0,2)

    self.params_svmWeights = QLineEdit(self.buttonGroup2,"params_svmWeights")
    self.params_svmWeights.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmWeights,2,2,0,1)

    self.params_svmDegree = QSpinBox(self.buttonGroup2,"params_svmDegree")
    self.params_svmDegree.setEnabled(0)
    self.params_svmDegree.setMinValue(1)
    self.params_svmDegree.setValue(3)

    layout22.addWidget(self.params_svmDegree,3,0)

    self.textLabel3 = QLabel(self.buttonGroup2,"textLabel3")

    layout22.addWidget(self.textLabel3,2,2)

    self.textLabel4_2 = QLabel(self.buttonGroup2,"textLabel4_2")

    layout22.addWidget(self.textLabel4_2,5,2)

    self.params_svmNu = QLineEdit(self.buttonGroup2,"params_svmNu")
    self.params_svmNu.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmNu,5,5,0,1)

    self.params_svmKernel = QComboBox(0,self.buttonGroup2,"params_svmKernel")
    self.params_svmKernel.setEnabled(0)

    layout22.addMultiCellWidget(self.params_svmKernel,0,0,0,1)
    layout23.addLayout(layout22)
    layout92.addLayout(layout23)
    layout93.addLayout(layout92)
    spacer60 = QSpacerItem(20,40,QSizePolicy.Minimum,QSizePolicy.Expanding)
    layout93.addItem(spacer60)
    buttonGroup2Layout.addLayout(layout93)
    TabPageLayout.addWidget(self.buttonGroup2)
    self.mainTab.insertTab(self.TabPage,QString(""))

    self.dbPage = QWidget(self.mainTab,"dbPage")
    self.mainTab.insertTab(self.dbPage,QString(""))

    self.descriptorPage = QWidget(self.mainTab,"descriptorPage")
    descriptorPageLayout = QVBoxLayout(self.descriptorPage,2,2,"descriptorPageLayout")

    self.descs_table = QTable(self.descriptorPage,"descs_table")
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.__tr("Name"))
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.__tr("Include"))
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.__tr("Role"))
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.__tr("Num Bounds"))
    self.descs_table.setNumRows(3)
    self.descs_table.setNumCols(4)
    descriptorPageLayout.addWidget(self.descs_table)

    Layout92 = QHBoxLayout(None,0,6,"Layout92")

    self.descs_refreshButton = QPushButton(self.descriptorPage,"descs_refreshButton")
    Layout92.addWidget(self.descs_refreshButton)
    Spacer25_2_2 = QSpacerItem(235,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout92.addItem(Spacer25_2_2)
    descriptorPageLayout.addLayout(Layout92)
    self.mainTab.insertTab(self.descriptorPage,QString(""))
    CompositeBuilderLayout.addWidget(self.mainTab)

    Layout13 = QHBoxLayout(None,0,6,"Layout13")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setAutoDefault(1)
    Layout13.addWidget(self.buttonHelp)
    Spacer14 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout13.addItem(Spacer14)

    self.buttonBuild = QPushButton(self,"buttonBuild")
    self.buttonBuild.setEnabled(0)
    Layout13.addWidget(self.buttonBuild)
    Spacer15 = QSpacerItem(16,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout13.addItem(Spacer15)

    self.buttonScreen = QPushButton(self,"buttonScreen")
    self.buttonScreen.setEnabled(0)
    Layout13.addWidget(self.buttonScreen)
    Spacer16 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout13.addItem(Spacer16)

    self.buttonInspect = QPushButton(self,"buttonInspect")
    self.buttonInspect.setEnabled(0)
    Layout13.addWidget(self.buttonInspect)
    Spacer15_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout13.addItem(Spacer15_2)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout13.addWidget(self.buttonClose)
    CompositeBuilderLayout.addLayout(Layout13)

    self.languageChange()

    self.resize(QSize(472,439).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.params_filter,SIGNAL("toggled(bool)"),self.toggleFilter)
    self.connect(self.params_split,SIGNAL("toggled(bool)"),self.toggleSplit)
    self.connect(self.buttonClose,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonBuild,SIGNAL("clicked()"),self.buildClick)
    self.connect(self.buttonScreen,SIGNAL("clicked()"),self.screenClick)
    self.connect(self.buttonHelp,SIGNAL("clicked()"),self.helpClick)
    self.connect(self.descs_refreshButton,SIGNAL("clicked()"),self.refreshClick)
    self.connect(self.params_lock,SIGNAL("toggled(bool)"),self.toggleLock)
    self.connect(self.buttonInspect,SIGNAL("clicked()"),self.inspectClick)
    self.connect(self.params_multiplebuilds,SIGNAL("toggled(bool)"),self.toggleMult)
    self.connect(self.params_treeRadio,SIGNAL("toggled(bool)"),self.selectTrees)
    self.connect(self.params_knnRadio,SIGNAL("toggled(bool)"),self.selectKNN)
    self.connect(self.params_svmRadio,SIGNAL("toggled(bool)"),self.selectSVM)
    self.connect(self.params_naiveBayesRadio,SIGNAL("toggled(bool)"),self.selectNaiveBayes)

    self.setTabOrder(self.params_count,self.params_bayes)
    self.setTabOrder(self.params_bayes,self.params_split)
    self.setTabOrder(self.params_split,self.params_splitFrac)
    self.setTabOrder(self.params_splitFrac,self.params_dataNone)
    self.setTabOrder(self.params_dataNone,self.params_actBounds)
    self.setTabOrder(self.params_actBounds,self.params_threshold)
    self.setTabOrder(self.params_threshold,self.params_multiplebuilds)
    self.setTabOrder(self.params_multiplebuilds,self.params_multiplebuildcount)
    self.setTabOrder(self.params_multiplebuildcount,self.params_filter)
    self.setTabOrder(self.params_filter,self.params_filterVal)
    self.setTabOrder(self.params_filterVal,self.params_filterFrac)
    self.setTabOrder(self.params_filterFrac,self.params_lock)
    self.setTabOrder(self.params_lock,self.params_lockV1)
    self.setTabOrder(self.params_lockV1,self.params_lockV2)
    self.setTabOrder(self.params_lockV2,self.params_persisttable)
    self.setTabOrder(self.params_persisttable,self.params_note)
    self.setTabOrder(self.params_note,self.params_treeRadio)
    self.setTabOrder(self.params_treeRadio,self.params_depth)
    self.setTabOrder(self.params_depth,self.params_autoBounds)
    self.setTabOrder(self.params_autoBounds,self.params_greedy)
    self.setTabOrder(self.params_greedy,self.params_recycle)
    self.setTabOrder(self.params_recycle,self.params_knnKVal)
    self.setTabOrder(self.params_knnKVal,self.params_knnMetric)
    self.setTabOrder(self.params_knnMetric,self.params_svmKernel)
    self.setTabOrder(self.params_svmKernel,self.params_svmType)
    self.setTabOrder(self.params_svmType,self.params_svmWeights)
    self.setTabOrder(self.params_svmWeights,self.params_svmDegree)
    self.setTabOrder(self.params_svmDegree,self.params_svmEps)
    self.setTabOrder(self.params_svmEps,self.params_svmNu)
    self.setTabOrder(self.params_svmNu,self.params_svmCoeff)
    self.setTabOrder(self.params_svmCoeff,self.params_svmGamma)
    self.setTabOrder(self.params_svmGamma,self.params_svmCost)
    self.setTabOrder(self.params_svmCost,self.buttonHelp)
    self.setTabOrder(self.buttonHelp,self.buttonBuild)
    self.setTabOrder(self.buttonBuild,self.buttonScreen)
    self.setTabOrder(self.buttonScreen,self.buttonInspect)
    self.setTabOrder(self.buttonInspect,self.buttonClose)
    self.setTabOrder(self.buttonClose,self.mainTab)
    self.setTabOrder(self.mainTab,self.descs_table)
    self.setTabOrder(self.descs_table,self.descs_refreshButton)


  def languageChange(self):
    self.setCaption(self.__tr("CompositeBuilder"))
    QToolTip.add(self.params_count,self.__tr("Number of individual models to include in the composite"))
    self.TextLabel4.setText(self.__tr("Number of Models"))
    QToolTip.add(self.TextLabel4,self.__tr("Number of individual models to include in the composite"))
    self.params_bayes.setText(self.__tr("Bayes?"))
    QToolTip.add(self.params_bayes,self.__tr("Toggles use of the Bayes composite algorithm."))
    self.GroupBox1.setTitle(self.__tr("Split Data"))
    self.params_split.setText(QString.null)
    QToolTip.add(self.params_split,self.__tr("Toggles splitting the data into training and holdout sets while building the composite."))
    self.params_splitFrac.setText(self.__tr("0.7"))
    QToolTip.add(self.params_splitFrac,self.__tr("Fraction of the data to be used in the training set."))
    self.TextLabel4_2_2_2.setText(self.__tr("Fraction"))
    QToolTip.add(self.TextLabel4_2_2_2,self.__tr("Fraction of the data to be used in the training set."))
    self.ButtonGroup1.setTitle(self.__tr("Randomization"))
    self.params_dataNone.setText(self.__tr("None"))
    QToolTip.add(self.params_dataNone,self.__tr("Do not alter the data"))
    self.params_dataShuffle.setText(self.__tr("Shuffle"))
    QToolTip.add(self.params_dataShuffle,self.__tr("Shuffle the activity values"))
    self.params_dataRandomize.setText(self.__tr("Randomize"))
    QToolTip.add(self.params_dataRandomize,self.__tr("Randomize the activity values"))
    self.params_actBounds.setText(QString.null)
    QToolTip.add(self.params_actBounds,self.__tr("quantization bounds to be used for activities"))
    self.TextLabel4_2_2_2_2.setText(self.__tr("Activity Bounds"))
    QToolTip.add(self.TextLabel4_2_2_2_2,self.__tr("quantization bounds to be used for activities"))
    self.params_threshold.setText(self.__tr("0.0"))
    QToolTip.add(self.params_threshold,self.__tr("Minimum confidence to be considered when making predictions."))
    self.TextLabel4_2_2.setText(self.__tr("Screen Threshold"))
    QToolTip.add(self.TextLabel4_2_2,self.__tr("Minimum confidence to be considered when making predictions."))
    self.GroupBox4.setTitle(self.__tr("Multiple Builds"))
    self.params_multiplebuilds.setText(self.__tr("Enable?"))
    QToolTip.add(self.params_multiplebuilds,self.__tr("Enables multiple builds (for statistics)."))
    QToolTip.add(self.params_multiplebuildcount,self.__tr("Number of composites to be built for the statistics."))
    self.TextLabel4_2_3.setText(self.__tr("Count"))
    QToolTip.add(self.TextLabel4_2_3,self.__tr("Maximum depth to be used when building individual models."))
    self.GroupBox2.setTitle(self.__tr("Filter Data"))
    self.params_filter.setText(self.__tr("Enable?"))
    QToolTip.add(self.params_filter,self.__tr("Toggles data filtering."))
    QToolTip.add(self.params_filterVal,self.__tr("The target value for filtering."))
    self.params_filterModels.setText(self.__tr("Models?"))
    QToolTip.add(self.params_filterModels,self.__tr("If checked, data filtering will be done at the model level instead of the full-dataset level."))
    self.TextLabel4_2_2_2_2_2.setText(self.__tr("Fraction"))
    QToolTip.add(self.TextLabel4_2_2_2_2_2,self.__tr("Target fraction for filtering."))
    self.TextLabel1.setText(self.__tr("Value"))
    QToolTip.add(self.TextLabel1,self.__tr("The target value for filtering."))
    self.params_filterFrac.setText(self.__tr("0.0"))
    QToolTip.add(self.params_filterFrac,self.__tr("Target fraction for filtering."))
    self.GroupBox3.setTitle(self.__tr("Lock Random"))
    self.params_lock.setText(self.__tr("Enable?"))
    QToolTip.add(self.params_lock,self.__tr("Enables locking the random number generator (mainly for debugging help)"))
    self.params_lockV2.setText(self.__tr("42"))
    QToolTip.add(self.params_lockV2,self.__tr("second number for the RNG lock"))
    self.TextLabel1_2.setText(self.__tr("V1"))
    self.TextLabel1_2_2.setText(self.__tr("V2"))
    self.params_lockV1.setText(self.__tr("23"))
    QToolTip.add(self.params_lockV1,self.__tr("first value for the RNG lock"))
    self.TextLabel1_3.setText(self.__tr("Persist Table"))
    self.TextLabel1_3_2.setText(self.__tr("Note"))
    QToolTip.add(self.params_note,self.__tr("Helpful text attached to the models when they are stored in the db."))
    QToolTip.add(self.params_persisttable,self.__tr("Database table in which to store the model(s) after being built."))
    self.mainTab.changeTab(self.paramPage,self.__tr("Parameters"))
    self.buttonGroup2.setTitle(self.__tr("Details"))
    self.params_treeRadio.setText(self.__tr("Trees"))
    QToolTip.add(self.params_treeRadio,self.__tr("Grow Trees"))
    QToolTip.add(self.params_depth,self.__tr("Maximum depth to be used when building individual models."))
    self.TextLabel4_2.setText(self.__tr("Max Depth"))
    QToolTip.add(self.TextLabel4_2,self.__tr("Maximum depth to be used when building individual models."))
    self.params_autoBounds.setText(self.__tr("Auto Bounds?"))
    QToolTip.add(self.params_autoBounds,self.__tr("Toggles use of the auto-bounds algorithm.  If this is not checked, you will need to provide quantization bounds for each descriptor."))
    self.params_greedy.setText(self.__tr("Less Greedy?"))
    QToolTip.add(self.params_greedy,self.__tr("Toggles use of a less greedy model-building algorithm."))
    self.params_recycle.setText(self.__tr("Recycle?"))
    QToolTip.add(self.params_recycle,self.__tr("Toggles recycling of descriptors (using them more than once in a single branch of a tree)."))
    self.params_prune.setText(self.__tr("Prune?"))
    QToolTip.add(self.params_prune,self.__tr("Toggles pruning"))
    QToolTip.add(self.params_randomDescriptors,self.__tr("If nonzero, random forests will be grown. Sets the number of descriptors available at each level of the random trees."))
    self.TextLabel4_2_4.setText(self.__tr("Random Vars?"))
    QToolTip.add(self.TextLabel4_2_4,self.__tr("If nonzero, random forests will be grown. Sets the number of descriptors available at each level of the random trees."))
    self.params_knnRadio.setText(self.__tr("KNN"))
    QToolTip.add(self.params_knnRadio,self.__tr("Use a KNN model"))
    QToolTip.add(self.params_knnKVal,self.__tr("The number of neighbors to consider"))
    self.textLabel1.setText(self.__tr("K"))
    QToolTip.add(self.textLabel1,self.__tr("The number of neighbors to consider"))
    self.params_knnMetric.clear()
    self.params_knnMetric.insertItem(self.__tr("Tanimoto"))
    self.params_knnMetric.insertItem(self.__tr("Euclidean"))
    self.params_knnMetric.setCurrentItem(0)
    QToolTip.add(self.params_knnMetric,self.__tr("Select the distance metric."))
    self.textLabel1_2.setText(self.__tr("Metric"))
    QToolTip.add(self.textLabel1_2,self.__tr("Select the distance metric."))
    self.params_naiveBayesRadio.setText(self.__tr("Naive Bayes"))
    QToolTip.add(self.params_naiveBayesRadio,self.__tr("Use a naive Bayes classifier"))
    self.params_naiveBayesM.setText(self.__tr("-1.0"))
    QToolTip.add(self.params_naiveBayesM,self.__tr("Enter the equivalent sample size (m in the m-estimate of probability). Values less than zero are ignored."))
    self.textLabel4_3.setText(self.__tr("Equiv. Sample Size"))
    QToolTip.add(self.textLabel4_3,self.__tr("Enter the equivalent sample size (m in the m-estimate of probability). Values less than zero are ignored."))
    self.params_svmRadio.setText(self.__tr("SVMs"))
    QToolTip.add(self.params_svmRadio,self.__tr("Use SVMs"))
    self.params_svmEps.setText(self.__tr("0.001"))
    QToolTip.add(self.params_svmEps,self.__tr("Enter the convergence criterion (epsilon) for the SVM building."))
    self.params_svmCost.setText(QString.null)
    QToolTip.add(self.params_svmCost,self.__tr("Enter the SVM cost value. Default is to find this via a grid search."))
    self.textLabel1_3.setText(self.__tr("Degree"))
    self.textLabel4_2_2.setText(self.__tr("Coeff"))
    self.params_svmGamma.setText(QString.null)
    QToolTip.add(self.params_svmGamma,self.__tr("Enter the SVM Gamma value. Default is to find this via a grid search."))
    self.params_svmCoeff.setText(self.__tr("0.0"))
    QToolTip.add(self.params_svmCoeff,self.__tr("Enter the coefficient (Coeff0) for the SVM."))
    self.textLabel6.setText(self.__tr("Cost"))
    self.textLabel5.setText(self.__tr("Gamma"))
    self.params_svmType.clear()
    self.params_svmType.insertItem(self.__tr("c-SVC"))
    self.params_svmType.insertItem(self.__tr("nu-SVC"))
    self.params_svmType.insertItem(self.__tr("one-class"))
    self.params_svmType.insertItem(self.__tr("eps-SVR"))
    self.params_svmType.insertItem(self.__tr("nu-SVR"))
    QToolTip.add(self.params_svmType,self.__tr("Select the classifier type."))
    self.textLabel4.setText(self.__tr("Epsilon"))
    self.textLabel2_2.setText(self.__tr("Type"))
    self.textLabel2.setText(self.__tr("Kernel"))
    QToolTip.add(self.params_svmWeights,self.__tr("Enter the weights for each possible activity as a comma-delimited list."))
    QToolTip.add(self.params_svmDegree,self.__tr("Degree of the kernel"))
    self.textLabel3.setText(self.__tr("Weights"))
    self.textLabel4_2.setText(self.__tr("Nu"))
    self.params_svmNu.setText(self.__tr("0.5"))
    QToolTip.add(self.params_svmNu,self.__tr("Enter the Nu value for the SVM."))
    self.params_svmKernel.clear()
    self.params_svmKernel.insertItem(self.__tr("linear"))
    self.params_svmKernel.insertItem(self.__tr("linearFP"))
    self.params_svmKernel.insertItem(self.__tr("polynomial"))
    self.params_svmKernel.insertItem(self.__tr("radial"))
    self.params_svmKernel.insertItem(self.__tr("radialFP"))
    self.params_svmKernel.insertItem(self.__tr("sigmoid"))
    self.params_svmKernel.setCurrentItem(3)
    QToolTip.add(self.params_svmKernel,self.__tr("Select the kernel type."))
    self.mainTab.changeTab(self.TabPage,self.__tr("Model Details"))
    self.mainTab.changeTab(self.dbPage,self.__tr("Database"))
    self.descs_table.horizontalHeader().setLabel(0,self.__tr("Name"))
    self.descs_table.horizontalHeader().setLabel(1,self.__tr("Include"))
    self.descs_table.horizontalHeader().setLabel(2,self.__tr("Role"))
    self.descs_table.horizontalHeader().setLabel(3,self.__tr("Num Bounds"))
    self.descs_refreshButton.setText(self.__tr("Refresh"))
    self.mainTab.changeTab(self.descriptorPage,self.__tr("Descriptors"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonHelp.setAccel(self.__tr("F1"))
    self.buttonBuild.setText(self.__tr("Build"))
    QToolTip.add(self.buttonBuild,self.__tr("Build a composite model."))
    self.buttonScreen.setText(self.__tr("Screen"))
    QToolTip.add(self.buttonScreen,self.__tr("Screen the composite model."))
    self.buttonInspect.setText(self.__tr("Inspect"))
    QToolTip.add(self.buttonInspect,self.__tr("Inspect the composite model (in another window)."))
    self.buttonClose.setText(self.__tr("Close"))
    self.buttonClose.setAccel(QString.null)


  def dbFileClick(self):
    print "CompositeBuilderDialog.dbFileClick(): Not implemented yet"

  def fileClick(self):
    print "CompositeBuilderDialog.fileClick(): Not implemented yet"

  def refreshClick(self):
    print "CompositeBuilderDialog.refreshClick(): Not implemented yet"

  def toggleFilter(self):
    print "CompositeBuilderDialog.toggleFilter(): Not implemented yet"

  def toggleSplit(self):
    print "CompositeBuilderDialog.toggleSplit(): Not implemented yet"

  def buildClick(self):
    print "CompositeBuilderDialog.buildClick(): Not implemented yet"

  def screenClick(self):
    print "CompositeBuilderDialog.screenClick(): Not implemented yet"

  def helpClick(self):
    print "CompositeBuilderDialog.helpClick(): Not implemented yet"

  def toggleLock(self):
    print "CompositeBuilderDialog.toggleLock(): Not implemented yet"

  def inspectClick(self):
    print "CompositeBuilderDialog.inspectClick(): Not implemented yet"

  def toggleMult(self):
    print "CompositeBuilderDialog.toggleMult(): Not implemented yet"

  def selectTrees(self):
    print "CompositeBuilderDialog.selectTrees(): Not implemented yet"

  def selectKNN(self):
    print "CompositeBuilderDialog.selectKNN(): Not implemented yet"

  def selectSVM(self):
    print "CompositeBuilderDialog.selectSVM(): Not implemented yet"

  def selectNaiveBayes(self):
    print "CompositeBuilderDialog.selectNaiveBayes(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("CompositeBuilderDialog",s,c)
