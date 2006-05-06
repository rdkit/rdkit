# Form implementation generated from reading ui file 'forms/LoaderWidget.ui'
#
# Created: Wed Nov 27 21:46:43 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *

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

class LoaderWidget(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    image0 = QPixmap(image0_data)

    if name == None:
      self.setName("LoaderWidget")

    self.resize(618,352)
    self.setCaption(self.trUtf8("Loader"))
    self.setIcon(image0)
    self.setSizeGripEnabled(0)

    LoaderWidgetLayout = QVBoxLayout(self,2,4,"LoaderWidgetLayout")

    self.supplierTab = QTabWidget(self,"supplierTab")

    self.infoPage = QWidget(self.supplierTab,"infoPage")
    infoPageLayout = QVBoxLayout(self.infoPage,2,4,"infoPageLayout")

    Layout72 = QHBoxLayout(None,0,6,"Layout72")

    self.info_typeGroup = QButtonGroup(self.infoPage,"info_typeGroup")
    self.info_typeGroup.setTitle(self.trUtf8("Input Type"))
    self.info_typeGroup.setExclusive(1)

    self.info_fileButton = QRadioButton(self.info_typeGroup,"info_fileButton")
    self.info_fileButton.setGeometry(QRect(10,50,130,24))
    self.info_fileButton.setText(self.trUtf8("File"))
    self.info_fileButton.setChecked(0)
    QToolTip.add(self.info_fileButton,self.trUtf8("Take compounds from a file.  See File tab."))

    self.info_molButton = QRadioButton(self.info_typeGroup,"info_molButton")
    self.info_molButton.setGeometry(QRect(10,20,130,24))
    self.info_molButton.setText(self.trUtf8("Molecule"))
    self.info_molButton.setChecked(0)
    QToolTip.add(self.info_molButton,self.trUtf8("Input a single compound.  See Molecule tab."))

    self.info_dbButton = QRadioButton(self.info_typeGroup,"info_dbButton")
    self.info_dbButton.setGeometry(QRect(10,80,130,24))
    self.info_dbButton.setText(self.trUtf8("Database"))
    self.info_dbButton.setChecked(1)
    QToolTip.add(self.info_dbButton,self.trUtf8("Take compounds from a database.  See Database tab."))
    Layout72.addWidget(self.info_typeGroup)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout72.addItem(spacer)
    infoPageLayout.addLayout(Layout72)

    Layout71 = QHBoxLayout(None,0,6,"Layout71")

    self.PushButton12 = QPushButton(self.infoPage,"PushButton12")
    self.PushButton12.setText(self.trUtf8("Refresh"))
    QToolTip.add(self.PushButton12,self.trUtf8("Read the file or database to preview the compounds."))
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

    self.molPage = QWidget(self.supplierTab,"molPage")
    molPageLayout = QVBoxLayout(self.molPage,2,4,"molPageLayout")

    self.mol_ImgButton = QToolButton(self.molPage,"mol_ImgButton")
    self.mol_ImgButton.setText(self.trUtf8("Click To Edit Molecule"))
    QToolTip.add(self.mol_ImgButton,self.trUtf8("Shows current molecule, click to edit."))
    molPageLayout.addWidget(self.mol_ImgButton)

    Layout33 = QHBoxLayout(None,0,6,"Layout33")

    self.TextLabel2 = QLabel(self.molPage,"TextLabel2")
    self.TextLabel2.setText(self.trUtf8("SMILES"))
    Layout33.addWidget(self.TextLabel2)

    self.mol_smilesBox = QLineEdit(self.molPage,"mol_smilesBox")
    QToolTip.add(self.mol_smilesBox,self.trUtf8("Enter SMILES for the molecule or by clicking above"))
    Layout33.addWidget(self.mol_smilesBox)
    molPageLayout.addLayout(Layout33)
    self.supplierTab.insertTab(self.molPage,self.trUtf8("Molecule"))

    self.filePage = QWidget(self.supplierTab,"filePage")
    filePageLayout = QVBoxLayout(self.filePage,11,6,"filePageLayout")

    Layout77 = QHBoxLayout(None,0,6,"Layout77")

    self.file_nameLabel = QLabel(self.filePage,"file_nameLabel")
    self.file_nameLabel.setText(self.trUtf8("File Name"))
    Layout77.addWidget(self.file_nameLabel)

    Layout76 = QHBoxLayout(None,0,6,"Layout76")

    self.file_nameBox = QLineEdit(self.filePage,"file_nameBox")
    QToolTip.add(self.file_nameBox,self.trUtf8("name of the file to use"))
    Layout76.addWidget(self.file_nameBox)

    self.file_nameButton = QToolButton(self.filePage,"file_nameButton")
    self.file_nameButton.setText(self.trUtf8("..."))
    QToolTip.add(self.file_nameButton,self.trUtf8("launches a file choice dialog box"))
    Layout76.addWidget(self.file_nameButton)
    Layout77.addLayout(Layout76)
    filePageLayout.addLayout(Layout77)

    Layout78 = QHBoxLayout(None,0,6,"Layout78")

    self.file_formatGroup = QButtonGroup(self.filePage,"file_formatGroup")
    self.file_formatGroup.setTitle(self.trUtf8("File Format"))
    self.file_formatGroup.setExclusive(1)
    self.file_formatGroup.setColumnLayout(0,Qt.Vertical)
    self.file_formatGroup.layout().setSpacing(6)
    self.file_formatGroup.layout().setMargin(11)
    file_formatGroupLayout = QVBoxLayout(self.file_formatGroup.layout())
    file_formatGroupLayout.setAlignment(Qt.AlignTop)

    self.file_sdfButton = QRadioButton(self.file_formatGroup,"file_sdfButton")
    self.file_sdfButton.setText(self.trUtf8("SD File"))
    self.file_sdfButton.setChecked(1)
    QToolTip.add(self.file_sdfButton,self.trUtf8("an MDL file format"))
    file_formatGroupLayout.addWidget(self.file_sdfButton)

    self.file_smiButton = QRadioButton(self.file_formatGroup,"file_smiButton")
    self.file_smiButton.setText(self.trUtf8("Smiles Table"))
    QToolTip.add(self.file_smiButton,self.trUtf8("a text file with SMILES in the first column"))
    file_formatGroupLayout.addWidget(self.file_smiButton)

    self.file_tdtButton = QRadioButton(self.file_formatGroup,"file_tdtButton")
    self.file_tdtButton.setText(self.trUtf8("TDT File"))
    QToolTip.add(self.file_tdtButton,self.trUtf8("a Daylight file format"))
    file_formatGroupLayout.addWidget(self.file_tdtButton)

    self.file_txtButton = QRadioButton(self.file_formatGroup,"file_txtButton")
    self.file_txtButton.setText(self.trUtf8("Raw Text"))
    QToolTip.add(self.file_txtButton,self.trUtf8("raw text file, the first line should contain column headers"))
    file_formatGroupLayout.addWidget(self.file_txtButton)
    Layout78.addWidget(self.file_formatGroup)
    spacer_4 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout78.addItem(spacer_4)
    filePageLayout.addLayout(Layout78)
    self.supplierTab.insertTab(self.filePage,self.trUtf8("File"))

    self.queryPage = QWidget(self.supplierTab,"queryPage")
    self.supplierTab.insertTab(self.queryPage,self.trUtf8("Database"))

    self.contentsPage = QWidget(self.supplierTab,"contentsPage")
    contentsPageLayout = QVBoxLayout(self.contentsPage,11,6,"contentsPageLayout")
    self.supplierTab.insertTab(self.contentsPage,self.trUtf8("Contents"))
    LoaderWidgetLayout.addWidget(self.supplierTab)

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
    LoaderWidgetLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.file_nameButton,SIGNAL("clicked()"),self.fileClick)
    self.connect(self.PushButton12,SIGNAL("clicked()"),self.refreshContents)

  def dbFileClick(self):
    print "LoaderWidget.dbFileClick(): Not implemented yet"

  def fileClick(self):
    print "LoaderWidget.fileClick(): Not implemented yet"

  def refreshContents(self):
    print "LoaderWidget.refreshContents(): Not implemented yet"
