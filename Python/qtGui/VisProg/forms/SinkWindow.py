# Form implementation generated from reading ui file 'SinkWindow.ui'
#
# Created: Fri Oct 25 14:54:36 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class SinkWindow(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("SinkWindow")

    self.resize(606,337)
    self.setCaption(self.trUtf8("Sink"))
    self.setSizeGripEnabled(0)

    SinkWindowLayout = QVBoxLayout(self,11,6,"SinkWindowLayout")

    self.sinkTab = QTabWidget(self,"sinkTab")

    self.infoPage = QWidget(self.sinkTab,"infoPage")
    infoPageLayout = QVBoxLayout(self.infoPage,11,6,"infoPageLayout")

    Layout30 = QVBoxLayout(None,0,6,"Layout30")

    Layout4 = QHBoxLayout(None,0,6,"Layout4")

    self.info_typeGroup = QButtonGroup(self.infoPage,"info_typeGroup")
    self.info_typeGroup.setTitle(self.trUtf8("Output Type"))
    self.info_typeGroup.setExclusive(1)

    self.info_dbButton = QRadioButton(self.info_typeGroup,"info_dbButton")
    self.info_dbButton.setGeometry(QRect(10,60,130,24))
    self.info_dbButton.setText(self.trUtf8("Database"))

    self.info_fileButton = QRadioButton(self.info_typeGroup,"info_fileButton")
    self.info_fileButton.setGeometry(QRect(10,30,130,24))
    self.info_fileButton.setText(self.trUtf8("File"))
    self.info_fileButton.setChecked(1)
    Layout4.addWidget(self.info_typeGroup)

    self.appendGroup = QButtonGroup(self.infoPage,"appendGroup")
    self.appendGroup.setTitle(self.trUtf8("Disposition"))

    self.info_appendButton = QRadioButton(self.appendGroup,"info_appendButton")
    self.info_appendButton.setGeometry(QRect(10,30,121,24))
    self.info_appendButton.setText(self.trUtf8("Append"))
    self.info_appendButton.setChecked(1)

    self.info_overwriteButton = QRadioButton(self.appendGroup,"info_overwriteButton")
    self.info_overwriteButton.setGeometry(QRect(10,60,121,24))
    self.info_overwriteButton.setText(self.trUtf8("Overwrite"))
    Layout4.addWidget(self.appendGroup)
    spacer = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout4.addItem(spacer)
    Layout30.addLayout(Layout4)

    Layout29 = QHBoxLayout(None,0,6,"Layout29")

    self.info_countBox = QSpinBox(self.infoPage,"info_countBox")
    self.info_countBox.setEnabled(1)
    self.info_countBox.setMaxValue(9999)
    QToolTip.add(self.info_countBox,self.trUtf8("Maximum number of items to provide"))
    Layout29.addWidget(self.info_countBox)

    self.TextLabel4 = QLabel(self.infoPage,"TextLabel4")
    self.TextLabel4.setText(self.trUtf8("Max Output Count"))
    QToolTip.add(self.TextLabel4,self.trUtf8("Maximum number of items to provide"))
    Layout29.addWidget(self.TextLabel4)
    spacer_2 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout29.addItem(spacer_2)
    Layout30.addLayout(Layout29)
    infoPageLayout.addLayout(Layout30)
    self.sinkTab.insertTab(self.infoPage,self.trUtf8("Info"))

    self.filePage = QWidget(self.sinkTab,"filePage")
    filePageLayout = QVBoxLayout(self.filePage,11,6,"filePageLayout")

    Layout6 = QHBoxLayout(None,0,6,"Layout6")

    self.file_nameLabel = QLabel(self.filePage,"file_nameLabel")
    self.file_nameLabel.setText(self.trUtf8("File Name"))
    Layout6.addWidget(self.file_nameLabel)

    self.file_nameBox = QLineEdit(self.filePage,"file_nameBox")
    Layout6.addWidget(self.file_nameBox)

    self.file_nameButton = QToolButton(self.filePage,"file_nameButton")
    self.file_nameButton.setText(self.trUtf8("..."))
    Layout6.addWidget(self.file_nameButton)
    filePageLayout.addLayout(Layout6)

    Layout27 = QHBoxLayout(None,0,6,"Layout27")

    self.file_formatGroup = QButtonGroup(self.filePage,"file_formatGroup")
    self.file_formatGroup.setTitle(self.trUtf8("File Format"))
    self.file_formatGroup.setExclusive(1)

    self.file_tdtButton = QRadioButton(self.file_formatGroup,"file_tdtButton")
    self.file_tdtButton.setGeometry(QRect(10,80,121,24))
    self.file_tdtButton.setText(self.trUtf8("TDT File"))

    self.file_txtButton = QRadioButton(self.file_formatGroup,"file_txtButton")
    self.file_txtButton.setGeometry(QRect(10,110,121,24))
    self.file_txtButton.setText(self.trUtf8("Raw Text"))

    self.file_sdfButton = QRadioButton(self.file_formatGroup,"file_sdfButton")
    self.file_sdfButton.setGeometry(QRect(10,20,121,24))
    self.file_sdfButton.setText(self.trUtf8("SD File"))
    self.file_sdfButton.setChecked(1)

    self.file_smiButton = QRadioButton(self.file_formatGroup,"file_smiButton")
    self.file_smiButton.setGeometry(QRect(10,50,121,24))
    self.file_smiButton.setText(self.trUtf8("Smiles Table"))
    Layout27.addWidget(self.file_formatGroup)
    spacer_3 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout27.addItem(spacer_3)
    filePageLayout.addLayout(Layout27)
    self.sinkTab.insertTab(self.filePage,self.trUtf8("File"))

    self.queryPage = QWidget(self.sinkTab,"queryPage")
    queryPageLayout = QGridLayout(self.queryPage,1,1,11,6,"queryPageLayout")

    self.db_sqlLabel_2 = QLabel(self.queryPage,"db_sqlLabel_2")
    self.db_sqlLabel_2.setText(self.trUtf8("SQL Where"))
    self.db_sqlLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    queryPageLayout.addWidget(self.db_sqlLabel_2,5,0)

    self.db_sqlWhereBox = QLineEdit(self.queryPage,"db_sqlWhereBox")
    QToolTip.add(self.db_sqlWhereBox,self.trUtf8("where specification for the SQL query"))

    queryPageLayout.addWidget(self.db_sqlWhereBox,5,1)

    self.db_sqlWhatBox = QLineEdit(self.queryPage,"db_sqlWhatBox")
    self.db_sqlWhatBox.setText(self.trUtf8("*"))
    QToolTip.add(self.db_sqlWhatBox,self.trUtf8("which fields are to be taken from the table "))

    queryPageLayout.addWidget(self.db_sqlWhatBox,4,1)

    self.db_tableCombo = QComboBox(0,self.queryPage,"db_tableCombo")
    self.db_tableCombo.setEnabled(0)
    self.db_tableCombo.setEditable(1)
    self.db_tableCombo.setAutoCompletion(1)
    self.db_tableCombo.setDuplicatesEnabled(0)
    QToolTip.add(self.db_tableCombo,self.trUtf8("name of the database table"))

    queryPageLayout.addWidget(self.db_tableCombo,3,1)

    self.db_tabelLabel = QLabel(self.queryPage,"db_tabelLabel")
    self.db_tabelLabel.setText(self.trUtf8("Table"))
    self.db_tabelLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    queryPageLayout.addWidget(self.db_tabelLabel,3,0)

    self.db_sqlLabel_3_2 = QLabel(self.queryPage,"db_sqlLabel_3_2")
    self.db_sqlLabel_3_2.setText(self.trUtf8("Password"))
    self.db_sqlLabel_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    QToolTip.add(self.db_sqlLabel_3_2,self.trUtf8("password for DB access"))

    queryPageLayout.addWidget(self.db_sqlLabel_3_2,2,0)

    self.db_sqlLabel_3 = QLabel(self.queryPage,"db_sqlLabel_3")
    self.db_sqlLabel_3.setText(self.trUtf8("User"))
    self.db_sqlLabel_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    QToolTip.add(self.db_sqlLabel_3,self.trUtf8("username for DB access"))

    queryPageLayout.addWidget(self.db_sqlLabel_3,1,0)

    self.db_passwordField = QLineEdit(self.queryPage,"db_passwordField")
    self.db_passwordField.setText(self.trUtf8("masterkey"))
    self.db_passwordField.setEchoMode(QLineEdit.Password)
    QToolTip.add(self.db_passwordField,self.trUtf8("password for DB access"))

    queryPageLayout.addWidget(self.db_passwordField,2,1)

    self.db_userField = QLineEdit(self.queryPage,"db_userField")
    self.db_userField.setText(self.trUtf8("sysdba"))
    QToolTip.add(self.db_userField,self.trUtf8("username for DB access"))

    queryPageLayout.addWidget(self.db_userField,1,1)

    self.db_sqlLabel = QLabel(self.queryPage,"db_sqlLabel")
    self.db_sqlLabel.setText(self.trUtf8("SQL What"))
    self.db_sqlLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    queryPageLayout.addWidget(self.db_sqlLabel,4,0)

    self.db_tabelLabel_2 = QLabel(self.queryPage,"db_tabelLabel_2")
    self.db_tabelLabel_2.setText(self.trUtf8("File"))
    self.db_tabelLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    queryPageLayout.addWidget(self.db_tabelLabel_2,0,0)

    Layout13 = QHBoxLayout(None,0,6,"Layout13")

    self.db_nameBox = QLineEdit(self.queryPage,"db_nameBox")
    QToolTip.add(self.db_nameBox,self.trUtf8("file containing the database"))
    Layout13.addWidget(self.db_nameBox)

    self.db_nameButton = QToolButton(self.queryPage,"db_nameButton")
    self.db_nameButton.setText(self.trUtf8("..."))
    Layout13.addWidget(self.db_nameButton)

    queryPageLayout.addLayout(Layout13,0,1)
    self.sinkTab.insertTab(self.queryPage,self.trUtf8("Database"))
    SinkWindowLayout.addWidget(self.sinkTab)

    Layout3 = QHBoxLayout(None,0,6,"Layout3")

    self.nodeTextLabel = QLabel(self,"nodeTextLabel")
    self.nodeTextLabel.setText(self.trUtf8("Node Text"))
    QToolTip.add(self.nodeTextLabel,self.trUtf8("Text displayed on the canvas"))
    Layout3.addWidget(self.nodeTextLabel)

    self.nodeTextBox = QLineEdit(self,"nodeTextBox")
    self.nodeTextBox.setText(self.trUtf8("Sink"))
    QToolTip.add(self.nodeTextBox,self.trUtf8("Text displayed on the canvas"))
    Layout3.addWidget(self.nodeTextBox)
    SinkWindowLayout.addLayout(Layout3)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer_4 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
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
    SinkWindowLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.db_nameButton,SIGNAL("clicked()"),self.dbFileClick)
    self.connect(self.file_nameButton,SIGNAL("clicked()"),self.fileClick)

  def dbFileClick(self):
    print "SinkWindow.dbFileClick(): Not implemented yet"

  def fileClick(self):
    print "SinkWindow.fileClick(): Not implemented yet"

  def refreshContents(self):
    print "SinkWindow.refreshContents(): Not implemented yet"
