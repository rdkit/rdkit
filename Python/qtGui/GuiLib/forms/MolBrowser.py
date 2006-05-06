# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MolBrowser.ui'
#
# Created: Mon Jan 23 08:36:05 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
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

class MolBrowser(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("MolBrowser")

    self.setIcon(self.image0)

    MolBrowserLayout = QVBoxLayout(self,2,4,"MolBrowserLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")

    self.inputPage = QWidget(self.tabWidget,"inputPage")
    inputPageLayout = QVBoxLayout(self.inputPage,2,4,"inputPageLayout")

    self.inputSourceGroup = QButtonGroup(self.inputPage,"inputSourceGroup")
    self.inputSourceGroup.setColumnLayout(0,Qt.Vertical)
    self.inputSourceGroup.layout().setSpacing(4)
    self.inputSourceGroup.layout().setMargin(2)
    inputSourceGroupLayout = QGridLayout(self.inputSourceGroup.layout())
    inputSourceGroupLayout.setAlignment(Qt.AlignTop)

    Layout79 = QHBoxLayout(None,0,6,"Layout79")

    self.inputFilenameBox = QLineEdit(self.inputSourceGroup,"inputFilenameBox")
    Layout79.addWidget(self.inputFilenameBox)

    self.inputFilenameButton = QToolButton(self.inputSourceGroup,"inputFilenameButton")
    Layout79.addWidget(self.inputFilenameButton)

    inputSourceGroupLayout.addLayout(Layout79,1,1)

    self.inputDbFrame = QFrame(self.inputSourceGroup,"inputDbFrame")
    self.inputDbFrame.setFrameShape(QFrame.Box)
    self.inputDbFrame.setFrameShadow(QFrame.Raised)

    inputSourceGroupLayout.addWidget(self.inputDbFrame,0,1)

    self.inputDbButton = QRadioButton(self.inputSourceGroup,"inputDbButton")
    self.inputDbButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.inputDbButton.sizePolicy().hasHeightForWidth()))
    self.inputDbButton.setChecked(1)

    inputSourceGroupLayout.addWidget(self.inputDbButton,0,0)

    self.inputFileButton = QRadioButton(self.inputSourceGroup,"inputFileButton")
    self.inputFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.inputFileButton.sizePolicy().hasHeightForWidth()))

    inputSourceGroupLayout.addWidget(self.inputFileButton,1,0)
    inputPageLayout.addWidget(self.inputSourceGroup)

    layout7 = QHBoxLayout(None,0,4,"layout7")

    self.inputLoadButton = QPushButton(self.inputPage,"inputLoadButton")
    layout7.addWidget(self.inputLoadButton)
    spacer4 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout7.addItem(spacer4)
    inputPageLayout.addLayout(layout7)
    self.tabWidget.insertTab(self.inputPage,QString.fromLatin1(""))

    self.dataPage = QWidget(self.tabWidget,"dataPage")
    self.tabWidget.insertTab(self.dataPage,QString.fromLatin1(""))

    self.outputPage = QWidget(self.tabWidget,"outputPage")

    LayoutWidget = QWidget(self.outputPage,"layout7_2")
    LayoutWidget.setGeometry(QRect(0,222,479,26))
    layout7_2 = QHBoxLayout(LayoutWidget,2,4,"layout7_2")

    self.outputSaveButton = QPushButton(LayoutWidget,"outputSaveButton")
    layout7_2.addWidget(self.outputSaveButton)
    spacer4_2 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout7_2.addItem(spacer4_2)

    self.outputSourceGroup = QButtonGroup(self.outputPage,"outputSourceGroup")
    self.outputSourceGroup.setGeometry(QRect(0,0,479,218))
    self.outputSourceGroup.setColumnLayout(0,Qt.Vertical)
    self.outputSourceGroup.layout().setSpacing(4)
    self.outputSourceGroup.layout().setMargin(2)
    outputSourceGroupLayout = QGridLayout(self.outputSourceGroup.layout())
    outputSourceGroupLayout.setAlignment(Qt.AlignTop)

    Layout79_2 = QHBoxLayout(None,0,6,"Layout79_2")

    self.outputFilenameBox = QLineEdit(self.outputSourceGroup,"outputFilenameBox")
    Layout79_2.addWidget(self.outputFilenameBox)

    self.outputFilenameButton = QToolButton(self.outputSourceGroup,"outputFilenameButton")
    Layout79_2.addWidget(self.outputFilenameButton)

    outputSourceGroupLayout.addLayout(Layout79_2,1,1)

    self.outputDbFrame = QFrame(self.outputSourceGroup,"outputDbFrame")
    self.outputDbFrame.setFrameShape(QFrame.Box)
    self.outputDbFrame.setFrameShadow(QFrame.Raised)

    outputSourceGroupLayout.addWidget(self.outputDbFrame,0,1)

    self.outputDbButton = QRadioButton(self.outputSourceGroup,"outputDbButton")
    self.outputDbButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.outputDbButton.sizePolicy().hasHeightForWidth()))
    self.outputDbButton.setChecked(1)

    outputSourceGroupLayout.addWidget(self.outputDbButton,0,0)

    self.outputFileButton = QRadioButton(self.outputSourceGroup,"outputFileButton")
    self.outputFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.outputFileButton.sizePolicy().hasHeightForWidth()))

    outputSourceGroupLayout.addWidget(self.outputFileButton,1,0)
    self.tabWidget.insertTab(self.outputPage,QString.fromLatin1(""))
    MolBrowserLayout.addWidget(self.tabWidget)

    self.languageChange()

    self.resize(QSize(491,287).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.inputLoadButton,SIGNAL("clicked()"),self.inputLoadClicked)
    self.connect(self.inputFilenameButton,SIGNAL("clicked()"),self.inputFileButtonClicked)
    self.connect(self.outputFilenameButton,SIGNAL("clicked()"),self.outputFileButtonClicked)
    self.connect(self.outputSaveButton,SIGNAL("clicked()"),self.outputSaveClicked)


  def languageChange(self):
    self.setCaption(self.__tr("MolBrowser"))
    self.inputSourceGroup.setTitle(self.__tr("Source"))
    QToolTip.add(self.inputFilenameBox,self.__tr("file containing the database"))
    self.inputFilenameButton.setText(self.__tr("..."))
    QToolTip.add(self.inputFilenameButton,self.__tr("Browse to the filename."))
    self.inputDbButton.setText(self.__tr("Database"))
    QToolTip.add(self.inputDbButton,self.__tr("Read molecules from a database table."))
    self.inputFileButton.setText(self.__tr("File"))
    QToolTip.add(self.inputFileButton,self.__tr("Read molecules from an input file (e.g. an SD or SMILES file)"))
    self.inputLoadButton.setText(self.__tr("Load"))
    QToolTip.add(self.inputLoadButton,self.__tr("Load the molecules into the Data page"))
    self.tabWidget.changeTab(self.inputPage,self.__tr("Input"))
    self.tabWidget.changeTab(self.dataPage,self.__tr("Data"))
    self.outputSaveButton.setText(self.__tr("Save"))
    QToolTip.add(self.outputSaveButton,self.__tr("Saves the molecules on the data page."))
    self.outputSourceGroup.setTitle(self.__tr("Destination"))
    QToolTip.add(self.outputFilenameBox,self.__tr("file containing the database"))
    self.outputFilenameButton.setText(self.__tr("..."))
    QToolTip.add(self.outputFilenameButton,self.__tr("Browse to the filename."))
    self.outputDbButton.setText(self.__tr("Database"))
    QToolTip.add(self.outputDbButton,self.__tr("Store the molecules in a database table."))
    self.outputFileButton.setText(self.__tr("File"))
    QToolTip.add(self.outputFileButton,self.__tr("Store molecules in an file (e.g. an SD or SMILES file)"))
    self.tabWidget.changeTab(self.outputPage,self.__tr("Output"))


  def inputLoadClicked(self):
    print "MolBrowser.inputLoadClicked(): Not implemented yet"

  def inputFileButtonClicked(self):
    print "MolBrowser.inputFileButtonClicked(): Not implemented yet"

  def outputSaveClicked(self):
    print "MolBrowser.outputSaveClicked(): Not implemented yet"

  def outputFileButtonClicked(self):
    print "MolBrowser.outputFileButtonClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("MolBrowser",s,c)
