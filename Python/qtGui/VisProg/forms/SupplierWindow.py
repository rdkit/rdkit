# Form implementation generated from reading ui file 'forms/SupplierWindow.ui'
#
# Created: Tue Oct 29 11:21:40 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *
from qttable import QTable


class SupplierWindow(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("SupplierWindow")

    self.resize(618,352)
    self.setCaption(self.trUtf8("Supplier"))
    self.setSizeGripEnabled(0)

    SupplierWindowLayout = QVBoxLayout(self,11,6,"SupplierWindowLayout")

    self.supplierTab = QTabWidget(self,"supplierTab")

    self.infoPage = QWidget(self.supplierTab,"infoPage")
    infoPageLayout = QVBoxLayout(self.infoPage,11,6,"infoPageLayout")

    Layout72 = QHBoxLayout(None,0,6,"Layout72")

    self.info_typeGroup = QButtonGroup(self.infoPage,"info_typeGroup")
    self.info_typeGroup.setTitle(self.trUtf8("Input Type"))
    self.info_typeGroup.setExclusive(1)

    self.info_dbButton = QRadioButton(self.info_typeGroup,"info_dbButton")
    self.info_dbButton.setGeometry(QRect(10,50,130,24))
    self.info_dbButton.setText(self.trUtf8("Database"))

    self.info_fileButton = QRadioButton(self.info_typeGroup,"info_fileButton")
    self.info_fileButton.setGeometry(QRect(10,20,130,24))
    self.info_fileButton.setText(self.trUtf8("File"))
    self.info_fileButton.setChecked(1)
    Layout72.addWidget(self.info_typeGroup)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout72.addItem(spacer)
    infoPageLayout.addLayout(Layout72)

    Layout71 = QHBoxLayout(None,0,6,"Layout71")

    self.PushButton12 = QPushButton(self.infoPage,"PushButton12")
    self.PushButton12.setText(self.trUtf8("Refresh"))
    Layout71.addWidget(self.PushButton12)
    spacer_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout71.addItem(spacer_2)

    self.info_countBox = QSpinBox(self.infoPage,"info_countBox")
    self.info_countBox.setEnabled(0)
    self.info_countBox.setMaxValue(99)
    QToolTip.add(self.info_countBox,self.trUtf8("Maximum number of items to provide"))
    Layout71.addWidget(self.info_countBox)

    self.TextLabel4 = QLabel(self.infoPage,"TextLabel4")
    self.TextLabel4.setText(self.trUtf8("Count"))
    QToolTip.add(self.TextLabel4,self.trUtf8("Maximum number of items to provide"))
    Layout71.addWidget(self.TextLabel4)
    spacer_3 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout71.addItem(spacer_3)
    infoPageLayout.addLayout(Layout71)
    self.supplierTab.insertTab(self.infoPage,self.trUtf8("Info"))

    self.filePage = QWidget(self.supplierTab,"filePage")
    filePageLayout = QVBoxLayout(self.filePage,11,6,"filePageLayout")

    Layout77 = QHBoxLayout(None,0,6,"Layout77")

    self.file_nameLabel = QLabel(self.filePage,"file_nameLabel")
    self.file_nameLabel.setText(self.trUtf8("File Name"))
    Layout77.addWidget(self.file_nameLabel)

    Layout76 = QHBoxLayout(None,0,6,"Layout76")

    self.file_nameBox = QLineEdit(self.filePage,"file_nameBox")
    Layout76.addWidget(self.file_nameBox)

    self.file_nameButton = QToolButton(self.filePage,"file_nameButton")
    self.file_nameButton.setText(self.trUtf8("..."))
    Layout76.addWidget(self.file_nameButton)
    Layout77.addLayout(Layout76)
    filePageLayout.addLayout(Layout77)

    Layout78 = QHBoxLayout(None,0,6,"Layout78")

    self.file_formatGroup = QButtonGroup(self.filePage,"file_formatGroup")
    self.file_formatGroup.setTitle(self.trUtf8("File Format"))
    self.file_formatGroup.setExclusive(1)

    self.file_smiButton = QRadioButton(self.file_formatGroup,"file_smiButton")
    self.file_smiButton.setGeometry(QRect(10,50,121,24))
    self.file_smiButton.setText(self.trUtf8("Smiles Table"))

    self.file_tdtButton = QRadioButton(self.file_formatGroup,"file_tdtButton")
    self.file_tdtButton.setGeometry(QRect(10,80,121,24))
    self.file_tdtButton.setText(self.trUtf8("TDT File"))

    self.file_sdfButton = QRadioButton(self.file_formatGroup,"file_sdfButton")
    self.file_sdfButton.setGeometry(QRect(10,20,121,24))
    self.file_sdfButton.setText(self.trUtf8("SD File"))
    self.file_sdfButton.setChecked(1)

    self.file_txtButton = QRadioButton(self.file_formatGroup,"file_txtButton")
    self.file_txtButton.setGeometry(QRect(10,110,121,24))
    self.file_txtButton.setText(self.trUtf8("Raw Text"))
    Layout78.addWidget(self.file_formatGroup)
    spacer_4 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout78.addItem(spacer_4)
    filePageLayout.addLayout(Layout78)
    self.supplierTab.insertTab(self.filePage,self.trUtf8("File"))

    self.queryPage = QWidget(self.supplierTab,"queryPage")
    self.supplierTab.insertTab(self.queryPage,self.trUtf8("Database"))

    self.contentsPage = QWidget(self.supplierTab,"contentsPage")
    contentsPageLayout = QVBoxLayout(self.contentsPage,11,6,"contentsPageLayout")

    self.supplierTable = QTable(self.contentsPage,"supplierTable")
    self.supplierTable.setNumRows(3)
    self.supplierTable.setNumCols(3)
    contentsPageLayout.addWidget(self.supplierTable)
    self.supplierTab.insertTab(self.contentsPage,self.trUtf8("Contents"))
    SupplierWindowLayout.addWidget(self.supplierTab)

    Layout74 = QHBoxLayout(None,0,6,"Layout74")

    self.nodeTextLabel = QLabel(self,"nodeTextLabel")
    self.nodeTextLabel.setText(self.trUtf8("Node Text"))
    QToolTip.add(self.nodeTextLabel,self.trUtf8("Text displayed on the canvas"))
    Layout74.addWidget(self.nodeTextLabel)

    self.nodeTextBox = QLineEdit(self,"nodeTextBox")
    self.nodeTextBox.setText(self.trUtf8("Supply"))
    QToolTip.add(self.nodeTextBox,self.trUtf8("Text displayed on the canvas"))
    Layout74.addWidget(self.nodeTextBox)
    SupplierWindowLayout.addLayout(Layout74)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer_5 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
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
    SupplierWindowLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.file_nameButton,SIGNAL("clicked()"),self.fileClick)
    self.connect(self.PushButton12,SIGNAL("clicked()"),self.refreshContents)

  def dbFileClick(self):
    print "SupplierWindow.dbFileClick(): Not implemented yet"

  def fileClick(self):
    print "SupplierWindow.fileClick(): Not implemented yet"

  def refreshContents(self):
    print "SupplierWindow.refreshContents(): Not implemented yet"
