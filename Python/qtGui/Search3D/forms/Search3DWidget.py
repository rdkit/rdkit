# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Search3DWidget.ui'
#
# Created: Thu Oct 19 08:53:50 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.16
#
# WARNING! All changes made in this file will be lost!


from qt import *


class Search3DWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if not name:
      self.setName("Search3DWidget")


    Search3DWidgetLayout = QVBoxLayout(self,4,2,"Search3DWidgetLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")

    self.loadPage = QWidget(self.tabWidget,"loadPage")
    loadPageLayout = QVBoxLayout(self.loadPage,4,2,"loadPageLayout")

    layout11 = QGridLayout(None,1,1,0,2,"layout11")

    self.proteinFilename = QLineEdit(self.loadPage,"proteinFilename")
    self.proteinFilename.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,QSizePolicy.Fixed,0,0,self.proteinFilename.sizePolicy().hasHeightForWidth()))
    self.proteinFilename.setMinimumSize(QSize(300,0))
    self.proteinFilename.setFrameShape(QLineEdit.LineEditPanel)
    self.proteinFilename.setFrameShadow(QLineEdit.Sunken)

    layout11.addWidget(self.proteinFilename,1,1)

    self.molFileButton = QPushButton(self.loadPage,"molFileButton")
    self.molFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.molFileButton.sizePolicy().hasHeightForWidth()))
    self.molFileButton.setMaximumSize(QSize(20,32767))

    layout11.addWidget(self.molFileButton,0,2)

    self.textLabel1_2 = QLabel(self.loadPage,"textLabel1_2")
    self.textLabel1_2.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.textLabel1_2.sizePolicy().hasHeightForWidth()))
    self.textLabel1_2.setMinimumSize(QSize(0,0))
    self.textLabel1_2.setScaledContents(0)
    self.textLabel1_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout11.addWidget(self.textLabel1_2,0,0)

    self.molFilename = QLineEdit(self.loadPage,"molFilename")
    self.molFilename.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,QSizePolicy.Fixed,0,0,self.molFilename.sizePolicy().hasHeightForWidth()))
    self.molFilename.setMinimumSize(QSize(300,0))
    self.molFilename.setFrameShape(QLineEdit.LineEditPanel)
    self.molFilename.setFrameShadow(QLineEdit.Sunken)

    layout11.addWidget(self.molFilename,0,1)

    self.proteinFileButton = QPushButton(self.loadPage,"proteinFileButton")
    self.proteinFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.proteinFileButton.sizePolicy().hasHeightForWidth()))
    self.proteinFileButton.setMaximumSize(QSize(20,32767))

    layout11.addWidget(self.proteinFileButton,1,2)

    self.textLabel1_2_2 = QLabel(self.loadPage,"textLabel1_2_2")
    self.textLabel1_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout11.addWidget(self.textLabel1_2_2,1,0)
    loadPageLayout.addLayout(layout11)
    spacer14 = QSpacerItem(20,131,QSizePolicy.Minimum,QSizePolicy.Expanding)
    loadPageLayout.addItem(spacer14)
    self.tabWidget.insertTab(self.loadPage,QString.fromLatin1(""))

    self.pcophorePage = QWidget(self.tabWidget,"pcophorePage")
    pcophorePageLayout = QVBoxLayout(self.pcophorePage,4,2,"pcophorePageLayout")

    layout13_2 = QHBoxLayout(None,0,2,"layout13_2")

    self.fdefFileLabel = QLabel(self.pcophorePage,"fdefFileLabel")
    self.fdefFileLabel.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.fdefFileLabel.sizePolicy().hasHeightForWidth()))
    self.fdefFileLabel.setMinimumSize(QSize(0,0))
    self.fdefFileLabel.setScaledContents(0)
    self.fdefFileLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    layout13_2.addWidget(self.fdefFileLabel)

    self.fdefFilename = QLineEdit(self.pcophorePage,"fdefFilename")
    self.fdefFilename.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,QSizePolicy.Fixed,0,0,self.fdefFilename.sizePolicy().hasHeightForWidth()))
    self.fdefFilename.setMinimumSize(QSize(300,0))
    self.fdefFilename.setMouseTracking(0)
    self.fdefFilename.setFrameShape(QLineEdit.LineEditPanel)
    self.fdefFilename.setFrameShadow(QLineEdit.Sunken)
    layout13_2.addWidget(self.fdefFilename)

    self.fdefFileButton = QPushButton(self.pcophorePage,"fdefFileButton")
    self.fdefFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.fdefFileButton.sizePolicy().hasHeightForWidth()))
    self.fdefFileButton.setMaximumSize(QSize(20,32767))
    layout13_2.addWidget(self.fdefFileButton)
    pcophorePageLayout.addLayout(layout13_2)

    layout7 = QHBoxLayout(None,0,2,"layout7")

    self.grabAtomsButton = QPushButton(self.pcophorePage,"grabAtomsButton")
    layout7.addWidget(self.grabAtomsButton)

    self.addAtomsButton = QPushButton(self.pcophorePage,"addAtomsButton")
    layout7.addWidget(self.addAtomsButton)
    pcophorePageLayout.addLayout(layout7)
    self.tabWidget.insertTab(self.pcophorePage,QString.fromLatin1(""))

    self.detailsPage = QWidget(self.tabWidget,"detailsPage")
    detailsPageLayout = QVBoxLayout(self.detailsPage,4,2,"detailsPageLayout")

    self.inputTypeGroup = QButtonGroup(self.detailsPage,"inputTypeGroup")
    self.inputTypeGroup.setColumnLayout(0,Qt.Vertical)
    self.inputTypeGroup.layout().setSpacing(2)
    self.inputTypeGroup.layout().setMargin(4)
    inputTypeGroupLayout = QVBoxLayout(self.inputTypeGroup.layout())
    inputTypeGroupLayout.setAlignment(Qt.AlignTop)

    layout9 = QHBoxLayout(None,0,2,"layout9")

    self.sdFileRadio = QRadioButton(self.inputTypeGroup,"sdFileRadio")
    self.sdFileRadio.setChecked(1)
    layout9.addWidget(self.sdFileRadio)

    self.molFileBox = QGroupBox(self.inputTypeGroup,"molFileBox")
    self.molFileBox.setColumnLayout(0,Qt.Vertical)
    self.molFileBox.layout().setSpacing(2)
    self.molFileBox.layout().setMargin(4)
    molFileBoxLayout = QHBoxLayout(self.molFileBox.layout())
    molFileBoxLayout.setAlignment(Qt.AlignTop)

    self.textLabel1_2_3 = QLabel(self.molFileBox,"textLabel1_2_3")
    self.textLabel1_2_3.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.textLabel1_2_3.sizePolicy().hasHeightForWidth()))
    self.textLabel1_2_3.setMinimumSize(QSize(0,0))
    self.textLabel1_2_3.setScaledContents(0)
    self.textLabel1_2_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    molFileBoxLayout.addWidget(self.textLabel1_2_3)

    self.sdFilename = QLineEdit(self.molFileBox,"sdFilename")
    self.sdFilename.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding,QSizePolicy.Fixed,0,0,self.sdFilename.sizePolicy().hasHeightForWidth()))
    self.sdFilename.setMinimumSize(QSize(300,0))
    self.sdFilename.setFrameShape(QLineEdit.LineEditPanel)
    self.sdFilename.setFrameShadow(QLineEdit.Sunken)
    molFileBoxLayout.addWidget(self.sdFilename)

    self.sdFileButton = QPushButton(self.molFileBox,"sdFileButton")
    self.sdFileButton.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.sdFileButton.sizePolicy().hasHeightForWidth()))
    self.sdFileButton.setMaximumSize(QSize(20,32767))
    molFileBoxLayout.addWidget(self.sdFileButton)
    layout9.addWidget(self.molFileBox)
    inputTypeGroupLayout.addLayout(layout9)
    detailsPageLayout.addWidget(self.inputTypeGroup)

    layout8 = QHBoxLayout(None,0,2,"layout8")

    self.loadMolsButton = QPushButton(self.detailsPage,"loadMolsButton")
    layout8.addWidget(self.loadMolsButton)

    self.addMolsButton = QPushButton(self.detailsPage,"addMolsButton")
    layout8.addWidget(self.addMolsButton)
    spacer4 = QSpacerItem(380,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout8.addItem(spacer4)
    detailsPageLayout.addLayout(layout8)

    self.numMolsLabel = QLabel(self.detailsPage,"numMolsLabel")
    detailsPageLayout.addWidget(self.numMolsLabel)

    self.searchButton = QPushButton(self.detailsPage,"searchButton")
    detailsPageLayout.addWidget(self.searchButton)

    self.numHitsLabel = QLabel(self.detailsPage,"numHitsLabel")
    detailsPageLayout.addWidget(self.numHitsLabel)

    self.examineHitsButton = QPushButton(self.detailsPage,"examineHitsButton")
    detailsPageLayout.addWidget(self.examineHitsButton)
    spacer7 = QSpacerItem(20,40,QSizePolicy.Minimum,QSizePolicy.Expanding)
    detailsPageLayout.addItem(spacer7)

    self.refineSearchButton = QPushButton(self.detailsPage,"refineSearchButton")
    detailsPageLayout.addWidget(self.refineSearchButton)

    layout7_2 = QHBoxLayout(None,0,2,"layout7_2")

    self.clusterHitsButton = QPushButton(self.detailsPage,"clusterHitsButton")
    layout7_2.addWidget(self.clusterHitsButton)

    self.textLabel1 = QLabel(self.detailsPage,"textLabel1")
    self.textLabel1.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.textLabel1.sizePolicy().hasHeightForWidth()))
    layout7_2.addWidget(self.textLabel1)

    self.numPicksBox = QSpinBox(self.detailsPage,"numPicksBox")
    self.numPicksBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.numPicksBox.sizePolicy().hasHeightForWidth()))
    layout7_2.addWidget(self.numPicksBox)
    detailsPageLayout.addLayout(layout7_2)
    self.tabWidget.insertTab(self.detailsPage,QString.fromLatin1(""))

    self.resultsPage = QWidget(self.tabWidget,"resultsPage")
    resultsPageLayout = QVBoxLayout(self.resultsPage,4,2,"resultsPageLayout")

    layout12 = QHBoxLayout(None,0,2,"layout12")

    self.centerViewCheck = QCheckBox(self.resultsPage,"centerViewCheck")
    self.centerViewCheck.setChecked(0)
    layout12.addWidget(self.centerViewCheck)

    self.replaceCheck = QCheckBox(self.resultsPage,"replaceCheck")
    self.replaceCheck.setChecked(1)
    layout12.addWidget(self.replaceCheck)

    self.showProteinCheck = QCheckBox(self.resultsPage,"showProteinCheck")
    self.showProteinCheck.setEnabled(1)
    layout12.addWidget(self.showProteinCheck)

    self.neighboringSurfaceCheck = QCheckBox(self.resultsPage,"neighboringSurfaceCheck")
    self.neighboringSurfaceCheck.setEnabled(1)
    layout12.addWidget(self.neighboringSurfaceCheck)

    self.showHBondsCheck = QCheckBox(self.resultsPage,"showHBondsCheck")
    layout12.addWidget(self.showHBondsCheck)

    self.showCollisionsCheck = QCheckBox(self.resultsPage,"showCollisionsCheck")
    layout12.addWidget(self.showCollisionsCheck)
    spacer9 = QSpacerItem(30,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout12.addItem(spacer9)
    resultsPageLayout.addLayout(layout12)
    self.tabWidget.insertTab(self.resultsPage,QString.fromLatin1(""))
    Search3DWidgetLayout.addWidget(self.tabWidget)

    self.languageChange()

    self.resize(QSize(647,510).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.sdFileButton,SIGNAL("clicked()"),self.sdFileButtonClicked)
    self.connect(self.sdFilename,SIGNAL("returnPressed()"),self.sdFilenameChanged)
    self.connect(self.searchButton,SIGNAL("clicked()"),self.searchButtonClicked)
    self.connect(self.molFilename,SIGNAL("returnPressed()"),self.molFilenameChanged)
    self.connect(self.molFileButton,SIGNAL("clicked()"),self.molFileButtonClicked)
    self.connect(self.proteinFilename,SIGNAL("returnPressed()"),self.proteinFilenameChanged)
    self.connect(self.proteinFileButton,SIGNAL("clicked()"),self.proteinFileButtonClicked)
    self.connect(self.fdefFilename,SIGNAL("returnPressed()"),self.fdefFilenameChanged)
    self.connect(self.fdefFileButton,SIGNAL("clicked()"),self.fdefFileButtonClicked)
    self.connect(self.grabAtomsButton,SIGNAL("clicked()"),self.grabAtomsButtonClicked)
    self.connect(self.loadMolsButton,SIGNAL("clicked()"),self.loadMolsButtonClicked)
    self.connect(self.examineHitsButton,SIGNAL("clicked()"),self.examineHitsClicked)
    self.connect(self.neighboringSurfaceCheck,SIGNAL("stateChanged(int)"),self.neighborStateChanged)
    self.connect(self.refineSearchButton,SIGNAL("clicked()"),self.refineSearchButtonClicked)
    self.connect(self.clusterHitsButton,SIGNAL("clicked()"),self.clusterHitsButtonClicked)
    self.connect(self.showProteinCheck,SIGNAL("stateChanged(int)"),self.showProteinStateChanged)
    self.connect(self.showHBondsCheck,SIGNAL("stateChanged(int)"),self.showHBondsStateChanged)
    self.connect(self.addAtomsButton,SIGNAL("clicked()"),self.addAtomsButtonClicked)
    self.connect(self.addMolsButton,SIGNAL("clicked()"),self.addMolsButtonClicked)
    self.connect(self.showCollisionsCheck,SIGNAL("stateChanged(int)"),self.showCollisionsStateChanged)


  def languageChange(self):
    self.setCaption(self.__tr("Search3DWidget"))
    self.proteinFilename.setInputMask(QString.null)
    QToolTip.add(self.proteinFilename,self.__tr("Enter the name of the mol file containing the reference PDB file. This is optional and is just used to display the protein."))
    self.molFileButton.setText(self.__tr("..."))
    QToolTip.add(self.molFileButton,self.__tr("Browse for the reference molecule file."))
    self.textLabel1_2.setText(self.__tr("Mol File:"))
    self.molFilename.setInputMask(QString.null)
    QToolTip.add(self.molFilename,self.__tr("Enter the name of the mol file containing the reference molecule."))
    self.proteinFileButton.setText(self.__tr("..."))
    QToolTip.add(self.proteinFileButton,self.__tr("Enter the name of the mol file containing the reference PDB file. This is optional and is just used to display the protein."))
    self.textLabel1_2_2.setText(self.__tr("PDB File (optional):"))
    self.tabWidget.changeTab(self.loadPage,self.__tr("Load Reference Molecule"))
    self.fdefFileLabel.setText(self.__tr("FDef File:"))
    self.fdefFilename.setInputMask(QString.null)
    QToolTip.add(self.fdefFilename,self.__tr("Enter the name of the feature definition file."))
    self.fdefFileButton.setText(self.__tr("..."))
    QToolTip.add(self.fdefFileButton,self.__tr("Browse for the feature definition file."))
    self.grabAtomsButton.setText(self.__tr("Grab Atoms From PyMol"))
    QToolTip.add(self.grabAtomsButton,self.__tr("Uses the current atom pick in PyMol to define a set of features."))
    self.addAtomsButton.setText(self.__tr("Add Atoms From PyMol"))
    QToolTip.add(self.addAtomsButton,self.__tr("Expands the feature set using the current selections in PyMol."))
    self.tabWidget.changeTab(self.pcophorePage,self.__tr("Define Pharmacophore"))
    self.inputTypeGroup.setTitle(QString.null)
    self.sdFileRadio.setText(QString.null)
    self.molFileBox.setTitle(QString.null)
    self.textLabel1_2_3.setText(self.__tr("SD File:"))
    self.sdFilename.setInputMask(QString.null)
    QToolTip.add(self.sdFilename,self.__tr("Enter the name of the mol file containing the SD file to be searched."))
    self.sdFileButton.setText(self.__tr("..."))
    QToolTip.add(self.sdFileButton,self.__tr("Browse for the SD file."))
    self.loadMolsButton.setText(self.__tr("Load Molecules"))
    QToolTip.add(self.loadMolsButton,self.__tr("Load molecules from the SD file. Replaces the current search set."))
    self.addMolsButton.setText(self.__tr("Add Molecules"))
    QToolTip.add(self.addMolsButton,self.__tr("Add molecules from the SD file to the current search set."))
    self.numMolsLabel.setText(self.__tr("Number of Molecules: 0"))
    self.searchButton.setText(self.__tr("Search"))
    QToolTip.add(self.searchButton,self.__tr("Search the target molecules to find those than can hit the active pharmacophore."))
    self.numHitsLabel.setText(self.__tr("Number of Hits:"))
    self.examineHitsButton.setText(self.__tr("Examine Hits"))
    QToolTip.add(self.examineHitsButton,self.__tr("Generate alignments for the search hits to the active pharmacophore."))
    self.refineSearchButton.setText(self.__tr("Refine Search"))
    QToolTip.add(self.refineSearchButton,self.__tr("Replaces the target molecules with the current search results (so that they can be searched again)."))
    self.clusterHitsButton.setText(self.__tr("Pick Diverse Subset"))
    QToolTip.add(self.clusterHitsButton,self.__tr("Cluster the hits and do a diversity pick."))
    self.textLabel1.setText(self.__tr("Number of picks:"))
    QToolTip.add(self.numPicksBox,self.__tr("Number of hits to include in the diverse subset."))
    self.tabWidget.changeTab(self.detailsPage,self.__tr("Search Details"))
    self.centerViewCheck.setText(self.__tr("Center View?"))
    QToolTip.add(self.centerViewCheck,self.__tr("If this is checked, the view will be recentered after each molecule is displayed."))
    self.replaceCheck.setText(self.__tr("Replace Current?"))
    QToolTip.add(self.replaceCheck,self.__tr("If this is checked, selecting a new alignment will replace the current alignment, otherwise both will be displayed."))
    self.showProteinCheck.setText(self.__tr("Show Protein?"))
    QToolTip.add(self.showProteinCheck,self.__tr("If this is checked and a protein file has been loaded, the protein will be displayed along with the ligand and new structures."))
    self.neighboringSurfaceCheck.setText(self.__tr("Show Neighborhood?"))
    QToolTip.add(self.neighboringSurfaceCheck,self.__tr("If this is checked, the surface of the protein in the neighborhood of the refrence molecule will be displayed."))
    self.showHBondsCheck.setText(self.__tr("Show H Bonds?"))
    QToolTip.add(self.showHBondsCheck,self.__tr("If this is checked, possible hydrogen bonds between the mapped molecule and the protein will be highlighted."))
    self.showCollisionsCheck.setText(self.__tr("Show Collisions?"))
    QToolTip.add(self.showCollisionsCheck,self.__tr("If this is checked, close contacts between the mapped molecule and the protein will be highlighted."))
    self.tabWidget.changeTab(self.resultsPage,self.__tr("Search Results"))


  def sdFileButtonClicked(self):
    print "Search3DWidget.sdFileButtonClicked(): Not implemented yet"

  def sdFilenameChanged(self):
    print "Search3DWidget.sdFilenameChanged(): Not implemented yet"

  def searchButtonClicked(self):
    print "Search3DWidget.searchButtonClicked(): Not implemented yet"

  def molFilenameChanged(self):
    print "Search3DWidget.molFilenameChanged(): Not implemented yet"

  def molFileButtonClicked(self):
    print "Search3DWidget.molFileButtonClicked(): Not implemented yet"

  def proteinFilenameChanged(self):
    print "Search3DWidget.proteinFilenameChanged(): Not implemented yet"

  def proteinFileButtonClicked(self):
    print "Search3DWidget.proteinFileButtonClicked(): Not implemented yet"

  def fdefFilenameChanged(self):
    print "Search3DWidget.fdefFilenameChanged(): Not implemented yet"

  def fdefFileButtonClicked(self):
    print "Search3DWidget.fdefFileButtonClicked(): Not implemented yet"

  def grabAtomsButtonClicked(self):
    print "Search3DWidget.grabAtomsButtonClicked(): Not implemented yet"

  def loadMolsButtonClicked(self):
    print "Search3DWidget.loadMolsButtonClicked(): Not implemented yet"

  def examineHitsClicked(self):
    print "Search3DWidget.examineHitsClicked(): Not implemented yet"

  def neighborStateChanged(self):
    print "Search3DWidget.neighborStateChanged(): Not implemented yet"

  def refineSearchButtonClicked(self):
    print "Search3DWidget.refineSearchButtonClicked(): Not implemented yet"

  def clusterHitsButtonClicked(self):
    print "Search3DWidget.clusterHitsButtonClicked(): Not implemented yet"

  def showProteinStateChanged(self):
    print "Search3DWidget.showProteinStateChanged(): Not implemented yet"

  def showHBondsStateChanged(self):
    print "Search3DWidget.showHBondsStateChanged(): Not implemented yet"

  def addAtomsButtonClicked(self):
    print "Search3DWidget.addAtomsButtonClicked(): Not implemented yet"

  def addMolsButtonClicked(self):
    print "Search3DWidget.addMolsButtonClicked(): Not implemented yet"

  def showCollisionsStateChanged(self):
    print "Search3DWidget.showCollisionsStateChanged(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("Search3DWidget",s,c)
