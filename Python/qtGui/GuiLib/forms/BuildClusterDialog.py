# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'BuildClusterDialog.ui'
#
# Created: Mon Jan 23 08:36:05 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
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

class ClusterBuilderDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("ClusterBuilderDialog")

    self.setIcon(self.image0)
    self.setSizeGripEnabled(0)

    ClusterBuilderDialogLayout = QVBoxLayout(self,2,4,"ClusterBuilderDialogLayout")

    self.mainTab = QTabWidget(self,"mainTab")

    self.paramPage = QWidget(self.mainTab,"paramPage")
    paramPageLayout = QVBoxLayout(self.paramPage,4,2,"paramPageLayout")

    layout21 = QHBoxLayout(None,0,2,"layout21")

    self.params_algorithmCombo = QComboBox(0,self.paramPage,"params_algorithmCombo")
    self.params_algorithmCombo.setMinimumSize(QSize(150,0))
    layout21.addWidget(self.params_algorithmCombo)

    self.TextLabel1_3 = QLabel(self.paramPage,"TextLabel1_3")
    layout21.addWidget(self.TextLabel1_3)
    spacer24 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout21.addItem(spacer24)
    paramPageLayout.addLayout(layout21)

    self.buttonGroup1 = QButtonGroup(self.paramPage,"buttonGroup1")
    self.buttonGroup1.setColumnLayout(0,Qt.Vertical)
    self.buttonGroup1.layout().setSpacing(2)
    self.buttonGroup1.layout().setMargin(4)
    buttonGroup1Layout = QGridLayout(self.buttonGroup1.layout())
    buttonGroup1Layout.setAlignment(Qt.AlignTop)

    self.paramsTypeDescriptors = QRadioButton(self.buttonGroup1,"paramsTypeDescriptors")
    self.paramsTypeDescriptors.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Fixed,0,0,self.paramsTypeDescriptors.sizePolicy().hasHeightForWidth()))
    self.paramsTypeDescriptors.setChecked(1)

    buttonGroup1Layout.addWidget(self.paramsTypeDescriptors,0,0)

    self.paramsTypeFingerprints = QRadioButton(self.buttonGroup1,"paramsTypeFingerprints")
    self.paramsTypeFingerprints.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Fixed,0,0,self.paramsTypeFingerprints.sizePolicy().hasHeightForWidth()))

    buttonGroup1Layout.addWidget(self.paramsTypeFingerprints,1,0)

    layout17 = QHBoxLayout(None,0,2,"layout17")

    self.params_fpMetricCombo = QComboBox(0,self.buttonGroup1,"params_fpMetricCombo")
    self.params_fpMetricCombo.setEnabled(0)
    layout17.addWidget(self.params_fpMetricCombo)

    self.TextLabel1_3_2_2 = QLabel(self.buttonGroup1,"TextLabel1_3_2_2")
    self.TextLabel1_3_2_2.setEnabled(1)
    self.TextLabel1_3_2_2.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Preferred,0,0,self.TextLabel1_3_2_2.sizePolicy().hasHeightForWidth()))
    layout17.addWidget(self.TextLabel1_3_2_2)

    buttonGroup1Layout.addLayout(layout17,1,1)

    layout15 = QGridLayout(None,1,1,0,2,"layout15")

    self.TextLabel1_3_2 = QLabel(self.buttonGroup1,"TextLabel1_3_2")
    self.TextLabel1_3_2.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Preferred,0,0,self.TextLabel1_3_2.sizePolicy().hasHeightForWidth()))
    self.TextLabel1_3_2.setFrameShape(QLabel.NoFrame)

    layout15.addWidget(self.TextLabel1_3_2,1,1)

    self.params_distanceCombo = QComboBox(0,self.buttonGroup1,"params_distanceCombo")
    self.params_distanceCombo.setMinimumSize(QSize(125,0))

    layout15.addWidget(self.params_distanceCombo,1,0)

    self.TextLabel1_3_3 = QLabel(self.buttonGroup1,"TextLabel1_3_3")
    self.TextLabel1_3_3.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Preferred,0,0,self.TextLabel1_3_3.sizePolicy().hasHeightForWidth()))

    layout15.addWidget(self.TextLabel1_3_3,0,1)

    self.params_standardCombo = QComboBox(0,self.buttonGroup1,"params_standardCombo")
    self.params_standardCombo.setMinimumSize(QSize(125,0))

    layout15.addWidget(self.params_standardCombo,0,0)

    buttonGroup1Layout.addLayout(layout15,0,1)
    spacer19 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    buttonGroup1Layout.addItem(spacer19,0,2)
    spacer19_2 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    buttonGroup1Layout.addItem(spacer19_2,1,2)
    paramPageLayout.addWidget(self.buttonGroup1)
    spacer25 = QSpacerItem(20,40,QSizePolicy.Minimum,QSizePolicy.Expanding)
    paramPageLayout.addItem(spacer25)

    layout19 = QHBoxLayout(None,0,2,"layout19")

    self.params_labelCheck = QCheckBox(self.paramPage,"params_labelCheck")
    self.params_labelCheck.setChecked(1)
    layout19.addWidget(self.params_labelCheck)
    Spacer12_2 = QSpacerItem(60,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout19.addItem(Spacer12_2)
    paramPageLayout.addLayout(layout19)

    layout20 = QHBoxLayout(None,0,2,"layout20")

    self.params_activityCheck = QCheckBox(self.paramPage,"params_activityCheck")
    self.params_activityCheck.setChecked(1)
    layout20.addWidget(self.params_activityCheck)
    Spacer12 = QSpacerItem(60,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout20.addItem(Spacer12)
    paramPageLayout.addLayout(layout20)

    layout5 = QHBoxLayout(None,0,2,"layout5")

    self.params_quantizeCheck = QCheckBox(self.paramPage,"params_quantizeCheck")
    layout5.addWidget(self.params_quantizeCheck)

    self.params_activityQuant = QLineEdit(self.paramPage,"params_activityQuant")
    self.params_activityQuant.setEnabled(0)
    self.params_activityQuant.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.params_activityQuant.sizePolicy().hasHeightForWidth()))
    layout5.addWidget(self.params_activityQuant)
    spacer7 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout5.addItem(spacer7)
    paramPageLayout.addLayout(layout5)
    self.mainTab.insertTab(self.paramPage,QString.fromLatin1(""))

    self.dbPage = QWidget(self.mainTab,"dbPage")
    self.mainTab.insertTab(self.dbPage,QString.fromLatin1(""))

    self.colsPage = QWidget(self.mainTab,"colsPage")

    LayoutWidget = QWidget(self.colsPage,"Layout92_2")
    LayoutWidget.setGeometry(QRect(10,280,573,32))
    Layout92_2 = QHBoxLayout(LayoutWidget,0,6,"Layout92_2")

    self.cols_refreshButton = QPushButton(LayoutWidget,"cols_refreshButton")
    Layout92_2.addWidget(self.cols_refreshButton)
    Spacer25_2_2_2 = QSpacerItem(235,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout92_2.addItem(Spacer25_2_2_2)

    self.cols_Table = QTable(self.colsPage,"cols_Table")
    self.cols_Table.setNumCols(self.cols_Table.numCols() + 1)
    self.cols_Table.horizontalHeader().setLabel(self.cols_Table.numCols() - 1,self.__tr("Name"))
    self.cols_Table.setNumCols(self.cols_Table.numCols() + 1)
    self.cols_Table.horizontalHeader().setLabel(self.cols_Table.numCols() - 1,self.__tr("Include?"))
    self.cols_Table.setNumCols(self.cols_Table.numCols() + 1)
    self.cols_Table.horizontalHeader().setLabel(self.cols_Table.numCols() - 1,self.__tr("Type"))
    self.cols_Table.setGeometry(QRect(10,10,560,260))
    self.cols_Table.setNumRows(0)
    self.cols_Table.setNumCols(3)
    self.mainTab.insertTab(self.colsPage,QString.fromLatin1(""))

    self.dataPage = QWidget(self.mainTab,"dataPage")
    self.mainTab.insertTab(self.dataPage,QString.fromLatin1(""))
    ClusterBuilderDialogLayout.addWidget(self.mainTab)

    Layout20 = QHBoxLayout(None,0,6,"Layout20")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setAutoDefault(1)
    Layout20.addWidget(self.buttonHelp)
    Spacer14 = QSpacerItem(40,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer14)

    self.refreshButton = QPushButton(self,"refreshButton")
    self.refreshButton.setEnabled(0)
    Layout20.addWidget(self.refreshButton)
    Spacer14_2 = QSpacerItem(40,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer14_2)

    self.buttonBuild = QPushButton(self,"buttonBuild")
    self.buttonBuild.setEnabled(0)
    Layout20.addWidget(self.buttonBuild)
    Spacer15_2 = QSpacerItem(70,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer15_2)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout20.addWidget(self.buttonClose)
    ClusterBuilderDialogLayout.addLayout(Layout20)

    self.languageChange()

    self.resize(QSize(432,286).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.buttonClose,SIGNAL("clicked()"),self.accept)
    self.connect(self.buttonBuild,SIGNAL("clicked()"),self.buildClick)
    self.connect(self.buttonHelp,SIGNAL("clicked()"),self.helpClick)
    self.connect(self.refreshButton,SIGNAL("clicked()"),self.refreshClick)
    self.connect(self.params_quantizeCheck,SIGNAL("stateChanged(int)"),self.quantizeChecked)
    self.connect(self.paramsTypeDescriptors,SIGNAL("toggled(bool)"),self.descriptorsChecked)
    self.connect(self.paramsTypeFingerprints,SIGNAL("toggled(bool)"),self.fingerprintsChecked)


  def languageChange(self):
    self.setCaption(self.__tr("ClusterBuilder"))
    QToolTip.add(self.params_algorithmCombo,self.__tr("Select the clustering algorithm to use."))
    self.TextLabel1_3.setText(self.__tr("Clustering Algorithm"))
    QToolTip.add(self.TextLabel1_3,self.__tr("Select the clustering algorithm to use."))
    self.buttonGroup1.setTitle(QString.null)
    self.paramsTypeDescriptors.setText(self.__tr("Descriptors"))
    QToolTip.add(self.paramsTypeDescriptors,self.__tr("The data are descriptor vectors."))
    self.paramsTypeFingerprints.setText(self.__tr("Fingerprints"))
    QToolTip.add(self.paramsTypeFingerprints,self.__tr("The data are fingerprints."))
    QToolTip.add(self.params_fpMetricCombo,self.__tr("Select the distance metric to use."))
    self.TextLabel1_3_2_2.setText(self.__tr("Distance Metric"))
    QToolTip.add(self.TextLabel1_3_2_2,self.__tr("Select the distance metric to use."))
    self.TextLabel1_3_2.setText(self.__tr("Distance Metric"))
    QToolTip.add(self.TextLabel1_3_2,self.__tr("Select the distance metric to use."))
    QToolTip.add(self.params_distanceCombo,self.__tr("Select the distance matrix to use."))
    self.TextLabel1_3_3.setText(self.__tr("Standardization"))
    QToolTip.add(self.TextLabel1_3_3,self.__tr("Select the standardization method to use."))
    QToolTip.add(self.params_standardCombo,self.__tr("Select the standardization to use."))
    self.params_labelCheck.setText(self.__tr("Labels Present?"))
    QToolTip.add(self.params_labelCheck,self.__tr("Set this toggle if there are activity labels in the data set to be clustered."))
    self.params_activityCheck.setText(self.__tr("Activities Present?"))
    QToolTip.add(self.params_activityCheck,self.__tr("Set this toggle if there are activity labels in the data set to be clustered."))
    self.params_quantizeCheck.setText(self.__tr("Quantize Activities?"))
    QToolTip.add(self.params_quantizeCheck,self.__tr("Use the provided quantization bound(s) to quantize the activities found in the db table."))
    QToolTip.add(self.params_activityQuant,self.__tr("Enter the quantization bound(s) for the activity"))
    self.mainTab.changeTab(self.paramPage,self.__tr("Parameters"))
    self.mainTab.changeTab(self.dbPage,self.__tr("Database"))
    self.cols_refreshButton.setText(self.__tr("Refresh"))
    QToolTip.add(self.cols_refreshButton,self.__tr("Refresh the column names from the database."))
    self.cols_Table.horizontalHeader().setLabel(0,self.__tr("Name"))
    self.cols_Table.horizontalHeader().setLabel(1,self.__tr("Include?"))
    self.cols_Table.horizontalHeader().setLabel(2,self.__tr("Type"))
    self.mainTab.changeTab(self.colsPage,self.__tr("Columns"))
    self.mainTab.changeTab(self.dataPage,self.__tr("Data"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonHelp.setAccel(self.__tr("F1"))
    self.refreshButton.setText(self.__tr("Refresh"))
    QToolTip.add(self.refreshButton,self.__tr("Refresh the data in the table from the database."))
    self.buttonBuild.setText(self.__tr("Build"))
    QToolTip.add(self.buttonBuild,self.__tr("Builds the cluster tree."))
    self.buttonClose.setText(self.__tr("Close"))
    self.buttonClose.setAccel(QString.null)


  def refreshClick(self):
    print "ClusterBuilderDialog.refreshClick(): Not implemented yet"

  def buildClick(self):
    print "ClusterBuilderDialog.buildClick(): Not implemented yet"

  def helpClick(self):
    print "ClusterBuilderDialog.helpClick(): Not implemented yet"

  def quantizeChecked(self):
    print "ClusterBuilderDialog.quantizeChecked(): Not implemented yet"

  def descriptorsChecked(self):
    print "ClusterBuilderDialog.descriptorsChecked(): Not implemented yet"

  def fingerprintsChecked(self):
    print "ClusterBuilderDialog.fingerprintsChecked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("ClusterBuilderDialog",s,c)
