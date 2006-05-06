# Form implementation generated from reading ui file '../forms/BuildCompositeDialog.ui'
#
# Created: Mon Oct 28 12:29:54 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *
from qttable import QTable


class CompositeBuilderDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("CompositeBuilder")

    self.resize(610,405)
    self.setCaption(self.trUtf8("CompositeBuilder"))
    self.setSizeGripEnabled(0)

    CompositeBuilderLayout = QVBoxLayout(self,11,6,"CompositeBuilderLayout")

    self.mainTab = QTabWidget(self,"mainTab")
    self.mainTab.setEnabled(1)

    self.paramPage = QWidget(self.mainTab,"paramPage")

    LayoutWidget = QWidget(self.paramPage,"Layout13")
    LayoutWidget.setGeometry(QRect(189,10,142,202))
    Layout13 = QVBoxLayout(LayoutWidget,0,6,"Layout13")

    self.GroupBox1 = QGroupBox(LayoutWidget,"GroupBox1")
    self.GroupBox1.setTitle(self.trUtf8("Split Data"))

    LayoutWidget_2 = QWidget(self.GroupBox1,"Layout19")
    LayoutWidget_2.setGeometry(QRect(11,24,118,25))
    Layout19 = QHBoxLayout(LayoutWidget_2,0,6,"Layout19")

    self.params_split = QCheckBox(LayoutWidget_2,"params_split")
    self.params_split.setText(self.trUtf8(""))
    self.params_split.setChecked(1)
    QToolTip.add(self.params_split,self.trUtf8("Toggles splitting the data into training and holdout sets while building the composite."))
    Layout19.addWidget(self.params_split)

    self.params_splitFrac = QLineEdit(LayoutWidget_2,"params_splitFrac")
    self.params_splitFrac.setSizePolicy(QSizePolicy(0,0,0,0,self.params_splitFrac.sizePolicy().hasHeightForWidth()))
    self.params_splitFrac.setMaximumSize(QSize(40,32767))
    self.params_splitFrac.setText(self.trUtf8("0.7"))
    QToolTip.add(self.params_splitFrac,self.trUtf8("Fraction of the data to be used in the training set."))
    Layout19.addWidget(self.params_splitFrac)

    self.TextLabel4_2_2_2 = QLabel(LayoutWidget_2,"TextLabel4_2_2_2")
    self.TextLabel4_2_2_2.setText(self.trUtf8("Fraction"))
    QToolTip.add(self.TextLabel4_2_2_2,self.trUtf8("Fraction of the data to be used in the training set."))
    Layout19.addWidget(self.TextLabel4_2_2_2)
    Layout13.addWidget(self.GroupBox1)

    self.ButtonGroup1 = QButtonGroup(LayoutWidget,"ButtonGroup1")
    self.ButtonGroup1.setTitle(self.trUtf8("Randomization"))

    self.params_dataShuffle = QRadioButton(self.ButtonGroup1,"params_dataShuffle")
    self.params_dataShuffle.setGeometry(QRect(10,50,64,21))
    self.params_dataShuffle.setText(self.trUtf8("Shuffle"))
    QToolTip.add(self.params_dataShuffle,self.trUtf8("Shuffle the activity values"))

    self.params_dataNone = QRadioButton(self.ButtonGroup1,"params_dataNone")
    self.params_dataNone.setGeometry(QRect(10,20,53,21))
    self.params_dataNone.setText(self.trUtf8("None"))
    self.params_dataNone.setChecked(1)
    QToolTip.add(self.params_dataNone,self.trUtf8("Do not alter the data"))

    self.params_dataRandomize = QRadioButton(self.ButtonGroup1,"params_dataRandomize")
    self.params_dataRandomize.setGeometry(QRect(10,80,101,21))
    self.params_dataRandomize.setText(self.trUtf8("Randomize"))
    QToolTip.add(self.params_dataRandomize,self.trUtf8("Randomize the activity values"))
    Layout13.addWidget(self.ButtonGroup1)

    LayoutWidget_3 = QWidget(self.paramPage,"Layout12")
    LayoutWidget_3.setGeometry(QRect(9,4,170,290))
    Layout12 = QGridLayout(LayoutWidget_3,1,1,0,6,"Layout12")

    self.params_greedy = QCheckBox(LayoutWidget_3,"params_greedy")
    self.params_greedy.setText(self.trUtf8("Less Greedy?"))
    QToolTip.add(self.params_greedy,self.trUtf8("Toggles use of a less greedy model-building algorithm."))

    Layout12.addMultiCellWidget(self.params_greedy,4,4,0,1)

    self.params_autoBounds = QCheckBox(LayoutWidget_3,"params_autoBounds")
    self.params_autoBounds.setText(self.trUtf8("Auto Bounds?"))
    self.params_autoBounds.setChecked(1)
    QToolTip.add(self.params_autoBounds,self.trUtf8("Toggles use of the auto-bounds algorithm.  If this is not checked, you will need to provide quantization bounds for each descriptor."))

    Layout12.addMultiCellWidget(self.params_autoBounds,3,3,0,1)

    self.TextLabel4_2 = QLabel(LayoutWidget_3,"TextLabel4_2")
    self.TextLabel4_2.setText(self.trUtf8("Max Depth"))
    QToolTip.add(self.TextLabel4_2,self.trUtf8("Maximum depth to be used when building individual models."))

    Layout12.addWidget(self.TextLabel4_2,1,1)

    self.TextLabel4 = QLabel(LayoutWidget_3,"TextLabel4")
    self.TextLabel4.setText(self.trUtf8("Number of Models"))
    QToolTip.add(self.TextLabel4,self.trUtf8("Number of individual models to include in the composite"))

    Layout12.addWidget(self.TextLabel4,0,1)

    self.params_threshold = QLineEdit(LayoutWidget_3,"params_threshold")
    self.params_threshold.setText(self.trUtf8("0.0"))
    QToolTip.add(self.params_threshold,self.trUtf8("Minimum confidence to be considered when making predictions."))

    Layout12.addWidget(self.params_threshold,7,0)
    spacer = QSpacerItem(0,57,QSizePolicy.Minimum,QSizePolicy.Expanding)
    Layout12.addItem(spacer,6,0)

    self.TextLabel4_2_2 = QLabel(LayoutWidget_3,"TextLabel4_2_2")
    self.TextLabel4_2_2.setText(self.trUtf8("Screen Threshold"))
    QToolTip.add(self.TextLabel4_2_2,self.trUtf8("Minimum confidence to be considered when making predictions."))

    Layout12.addWidget(self.TextLabel4_2_2,7,1)

    self.params_depth = QSpinBox(LayoutWidget_3,"params_depth")
    self.params_depth.setMaxValue(1000)
    self.params_depth.setMinValue(-1)
    self.params_depth.setValue(-1)
    QToolTip.add(self.params_depth,self.trUtf8("Maximum depth to be used when building individual models."))

    Layout12.addWidget(self.params_depth,1,0)

    self.params_count = QSpinBox(LayoutWidget_3,"params_count")
    self.params_count.setMaxValue(1000)
    self.params_count.setLineStep(5)
    self.params_count.setValue(10)
    QToolTip.add(self.params_count,self.trUtf8("Number of individual models to include in the composite"))

    Layout12.addWidget(self.params_count,0,0)
    spacer_2 = QSpacerItem(0,57,QSizePolicy.Minimum,QSizePolicy.Expanding)
    Layout12.addItem(spacer_2,2,0)

    self.params_bayes = QCheckBox(LayoutWidget_3,"params_bayes")
    self.params_bayes.setText(self.trUtf8("Bayes?"))
    self.params_bayes.setChecked(0)
    QToolTip.add(self.params_bayes,self.trUtf8("Toggles use of the Bayes composite algorithm."))

    Layout12.addMultiCellWidget(self.params_bayes,5,5,0,1)

    LayoutWidget_4 = QWidget(self.paramPage,"Layout22")
    LayoutWidget_4.setGeometry(QRect(420,10,167,254))
    Layout22 = QVBoxLayout(LayoutWidget_4,0,6,"Layout22")

    self.GroupBox3 = QGroupBox(LayoutWidget_4,"GroupBox3")
    self.GroupBox3.setTitle(self.trUtf8("Lock Random"))

    LayoutWidget_5 = QWidget(self.GroupBox3,"Layout21")
    LayoutWidget_5.setGeometry(QRect(8,24,132,83))
    Layout21 = QVBoxLayout(LayoutWidget_5,0,6,"Layout21")

    self.params_lock = QCheckBox(LayoutWidget_5,"params_lock")
    self.params_lock.setText(self.trUtf8("Enable?"))
    QToolTip.add(self.params_lock,self.trUtf8("Enables locking the random number generator (mainly for debugging help)"))
    Layout21.addWidget(self.params_lock)

    Layout20 = QGridLayout(None,1,1,0,6,"Layout20")

    self.params_lockV2 = QLineEdit(LayoutWidget_5,"params_lockV2")
    self.params_lockV2.setEnabled(0)
    self.params_lockV2.setMaximumSize(QSize(40,32767))
    self.params_lockV2.setText(self.trUtf8("42"))
    QToolTip.add(self.params_lockV2,self.trUtf8("second number for the RNG lock"))

    Layout20.addWidget(self.params_lockV2,1,1)
    spacer_3 = QSpacerItem(56,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(spacer_3,1,0)

    self.TextLabel1_2 = QLabel(LayoutWidget_5,"TextLabel1_2")
    self.TextLabel1_2.setText(self.trUtf8("V1"))

    Layout20.addWidget(self.TextLabel1_2,0,2)

    self.TextLabel1_2_2 = QLabel(LayoutWidget_5,"TextLabel1_2_2")
    self.TextLabel1_2_2.setText(self.trUtf8("V2"))

    Layout20.addWidget(self.TextLabel1_2_2,1,2)

    self.params_lockV1 = QLineEdit(LayoutWidget_5,"params_lockV1")
    self.params_lockV1.setEnabled(0)
    self.params_lockV1.setMaximumSize(QSize(40,32767))
    self.params_lockV1.setText(self.trUtf8("23"))
    QToolTip.add(self.params_lockV1,self.trUtf8("first value for the RNG lock"))

    Layout20.addWidget(self.params_lockV1,0,1)
    spacer_4 = QSpacerItem(56,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(spacer_4,0,0)
    Layout21.addLayout(Layout20)
    Layout22.addWidget(self.GroupBox3)

    self.GroupBox2 = QGroupBox(LayoutWidget_4,"GroupBox2")
    self.GroupBox2.setTitle(self.trUtf8("Filter Data"))

    LayoutWidget_6 = QWidget(self.GroupBox2,"Layout19")
    LayoutWidget_6.setGeometry(QRect(11,16,143,83))
    Layout19_2 = QVBoxLayout(LayoutWidget_6,0,6,"Layout19_2")

    self.params_filter = QCheckBox(LayoutWidget_6,"params_filter")
    self.params_filter.setText(self.trUtf8("Enable?"))
    QToolTip.add(self.params_filter,self.trUtf8("Toggles data filtering."))
    Layout19_2.addWidget(self.params_filter)

    Layout18 = QGridLayout(None,1,1,0,6,"Layout18")

    self.params_filterFrac = QLineEdit(LayoutWidget_6,"params_filterFrac")
    self.params_filterFrac.setEnabled(0)
    self.params_filterFrac.setSizePolicy(QSizePolicy(0,0,0,0,self.params_filterFrac.sizePolicy().hasHeightForWidth()))
    self.params_filterFrac.setMaximumSize(QSize(40,32767))
    self.params_filterFrac.setText(self.trUtf8("0.0"))
    QToolTip.add(self.params_filterFrac,self.trUtf8("Target fraction for filtering."))

    Layout18.addWidget(self.params_filterFrac,1,1)

    self.TextLabel4_2_2_2_2 = QLabel(LayoutWidget_6,"TextLabel4_2_2_2_2")
    self.TextLabel4_2_2_2_2.setText(self.trUtf8("Fraction"))
    QToolTip.add(self.TextLabel4_2_2_2_2,self.trUtf8("Target fraction for filtering."))

    Layout18.addWidget(self.TextLabel4_2_2_2_2,1,2)
    spacer_5 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout18.addItem(spacer_5,0,0)
    spacer_6 = QSpacerItem(16,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout18.addItem(spacer_6,1,0)

    self.params_filterVal = QSpinBox(LayoutWidget_6,"params_filterVal")
    self.params_filterVal.setEnabled(0)
    QToolTip.add(self.params_filterVal,self.trUtf8("The target value for filtering."))

    Layout18.addWidget(self.params_filterVal,0,1)

    self.TextLabel1 = QLabel(LayoutWidget_6,"TextLabel1")
    self.TextLabel1.setText(self.trUtf8("Value"))
    QToolTip.add(self.TextLabel1,self.trUtf8("The target value for filtering."))

    Layout18.addWidget(self.TextLabel1,0,2)
    Layout19_2.addLayout(Layout18)
    Layout22.addWidget(self.GroupBox2)
    self.mainTab.insertTab(self.paramPage,self.trUtf8("Parameters"))

    self.dbPage = QWidget(self.mainTab,"dbPage")
    dbPageLayout = QVBoxLayout(self.dbPage,11,6,"dbPageLayout")

    Layout13_2 = QGridLayout(None,1,1,0,6,"Layout13_2")

    self.db_sqlLabel_3 = QLabel(self.dbPage,"db_sqlLabel_3")
    self.db_sqlLabel_3.setText(self.trUtf8("User"))
    self.db_sqlLabel_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    QToolTip.add(self.db_sqlLabel_3,self.trUtf8("username for DB access"))

    Layout13_2.addWidget(self.db_sqlLabel_3,1,0)

    self.db_tableCombo = QComboBox(0,self.dbPage,"db_tableCombo")
    self.db_tableCombo.setEnabled(0)
    self.db_tableCombo.setAutoCompletion(1)
    self.db_tableCombo.setDuplicatesEnabled(0)
    QToolTip.add(self.db_tableCombo,self.trUtf8("name of the database table"))

    Layout13_2.addWidget(self.db_tableCombo,3,1)

    self.db_passwordField = QLineEdit(self.dbPage,"db_passwordField")
    self.db_passwordField.setText(self.trUtf8("masterkey"))
    self.db_passwordField.setEchoMode(QLineEdit.Password)
    QToolTip.add(self.db_passwordField,self.trUtf8("password for DB access"))

    Layout13_2.addWidget(self.db_passwordField,2,1)

    self.db_tabelLabel_2 = QLabel(self.dbPage,"db_tabelLabel_2")
    self.db_tabelLabel_2.setText(self.trUtf8("File"))
    self.db_tabelLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout13_2.addWidget(self.db_tabelLabel_2,0,0)

    self.db_tabelLabel = QLabel(self.dbPage,"db_tabelLabel")
    self.db_tabelLabel.setText(self.trUtf8("Table"))
    self.db_tabelLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout13_2.addWidget(self.db_tabelLabel,3,0)

    self.db_sqlJoinBox = QLineEdit(self.dbPage,"db_sqlJoinBox")
    QToolTip.add(self.db_sqlJoinBox,self.trUtf8("where specification for the SQL query"))

    Layout13_2.addWidget(self.db_sqlJoinBox,6,1)

    self.db_sqlWhatBox = QLineEdit(self.dbPage,"db_sqlWhatBox")
    self.db_sqlWhatBox.setText(self.trUtf8("*"))
    QToolTip.add(self.db_sqlWhatBox,self.trUtf8("which fields are to be taken from the table "))

    Layout13_2.addWidget(self.db_sqlWhatBox,4,1)

    self.db_sqlLabel_2_2 = QLabel(self.dbPage,"db_sqlLabel_2_2")
    self.db_sqlLabel_2_2.setText(self.trUtf8("SQL Join"))
    self.db_sqlLabel_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout13_2.addWidget(self.db_sqlLabel_2_2,6,0)

    self.db_sqlLabel_2 = QLabel(self.dbPage,"db_sqlLabel_2")
    self.db_sqlLabel_2.setText(self.trUtf8("SQL Where"))
    self.db_sqlLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout13_2.addWidget(self.db_sqlLabel_2,5,0)

    self.db_sqlLabel = QLabel(self.dbPage,"db_sqlLabel")
    self.db_sqlLabel.setText(self.trUtf8("SQL What"))
    self.db_sqlLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    Layout13_2.addWidget(self.db_sqlLabel,4,0)

    Layout79 = QHBoxLayout(None,0,6,"Layout79")

    self.db_nameBox = QLineEdit(self.dbPage,"db_nameBox")
    QToolTip.add(self.db_nameBox,self.trUtf8("file containing the database"))
    Layout79.addWidget(self.db_nameBox)

    self.db_nameButton = QToolButton(self.dbPage,"db_nameButton")
    self.db_nameButton.setText(self.trUtf8("..."))
    Layout79.addWidget(self.db_nameButton)

    Layout13_2.addLayout(Layout79,0,1)
    spacer_7 = QSpacerItem(0,20,QSizePolicy.Minimum,QSizePolicy.Expanding)
    Layout13_2.addItem(spacer_7,7,1)

    self.db_sqlLabel_3_2 = QLabel(self.dbPage,"db_sqlLabel_3_2")
    self.db_sqlLabel_3_2.setText(self.trUtf8("Password"))
    self.db_sqlLabel_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    QToolTip.add(self.db_sqlLabel_3_2,self.trUtf8("password for DB access"))

    Layout13_2.addWidget(self.db_sqlLabel_3_2,2,0)

    self.db_userField = QLineEdit(self.dbPage,"db_userField")
    self.db_userField.setText(self.trUtf8("sysdba"))
    QToolTip.add(self.db_userField,self.trUtf8("username for DB access"))

    Layout13_2.addWidget(self.db_userField,1,1)

    self.db_sqlWhereBox = QLineEdit(self.dbPage,"db_sqlWhereBox")
    QToolTip.add(self.db_sqlWhereBox,self.trUtf8("where specification for the SQL query"))

    Layout13_2.addWidget(self.db_sqlWhereBox,5,1)
    dbPageLayout.addLayout(Layout13_2)
    self.mainTab.insertTab(self.dbPage,self.trUtf8("Database"))

    self.descriptorPage = QWidget(self.mainTab,"descriptorPage")
    descriptorPageLayout = QVBoxLayout(self.descriptorPage,11,6,"descriptorPageLayout")

    self.descs_table = QTable(self.descriptorPage,"descs_table")
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.trUtf8("Name"))
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.trUtf8("Role"))
    self.descs_table.setNumCols(self.descs_table.numCols() + 1)
    self.descs_table.horizontalHeader().setLabel(self.descs_table.numCols() - 1,self.trUtf8("Num Bounds"))
    self.descs_table.setNumRows(3)
    self.descs_table.setNumCols(3)
    descriptorPageLayout.addWidget(self.descs_table)

    Layout92 = QHBoxLayout(None,0,6,"Layout92")

    self.descs_refreshButton = QPushButton(self.descriptorPage,"descs_refreshButton")
    self.descs_refreshButton.setText(self.trUtf8("Refresh"))
    Layout92.addWidget(self.descs_refreshButton)
    spacer_8 = QSpacerItem(235,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout92.addItem(spacer_8)
    descriptorPageLayout.addLayout(Layout92)
    self.mainTab.insertTab(self.descriptorPage,self.trUtf8("Descriptors"))
    CompositeBuilderLayout.addWidget(self.mainTab)

    Layout19_3 = QHBoxLayout(None,0,6,"Layout19_3")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout19_3.addWidget(self.buttonHelp)
    spacer_9 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout19_3.addItem(spacer_9)

    self.buttonBuild = QPushButton(self,"buttonBuild")
    self.buttonBuild.setEnabled(0)
    self.buttonBuild.setText(self.trUtf8("Build"))
    QToolTip.add(self.buttonBuild,self.trUtf8("Build a composite model."))
    Layout19_3.addWidget(self.buttonBuild)
    spacer_10 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout19_3.addItem(spacer_10)

    self.buttonScreen = QPushButton(self,"buttonScreen")
    self.buttonScreen.setEnabled(0)
    self.buttonScreen.setText(self.trUtf8("Screen"))
    QToolTip.add(self.buttonScreen,self.trUtf8("Build a composite model."))
    Layout19_3.addWidget(self.buttonScreen)
    spacer_11 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout19_3.addItem(spacer_11)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setText(self.trUtf8("Close"))
    self.buttonClose.setAccel(0)
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout19_3.addWidget(self.buttonClose)
    CompositeBuilderLayout.addLayout(Layout19_3)

    self.connect(self.params_filter,SIGNAL("toggled(bool)"),self.toggleFilter)
    self.connect(self.params_split,SIGNAL("toggled(bool)"),self.toggleSplit)
    self.connect(self.buttonClose,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonBuild,SIGNAL("clicked()"),self.buildClick)
    self.connect(self.buttonScreen,SIGNAL("clicked()"),self.screenClick)
    self.connect(self.buttonHelp,SIGNAL("clicked()"),self.helpClick)
    self.connect(self.descs_refreshButton,SIGNAL("clicked()"),self.refreshClick)
    self.connect(self.db_nameButton,SIGNAL("clicked()"),self.dbFileClick)
    self.connect(self.params_lock,SIGNAL("toggled(bool)"),self.toggleLock)

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
