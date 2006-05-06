# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MolDescriptorWin.ui'
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
image1_data = [
"22 22 5 1",
". c None",
"# c #000000",
"c c #848200",
"a c #ffff00",
"b c #ffffff",
"......................",
"......................",
"......................",
"............####....#.",
"...........#....##.##.",
"..................###.",
".................####.",
".####...........#####.",
"#abab##########.......",
"#babababababab#.......",
"#ababababababa#.......",
"#babababababab#.......",
"#ababab###############",
"#babab##cccccccccccc##",
"#abab##cccccccccccc##.",
"#bab##cccccccccccc##..",
"#ab##cccccccccccc##...",
"#b##cccccccccccc##....",
"###cccccccccccc##.....",
"##cccccccccccc##......",
"###############.......",
"......................"
]
image2_data = [
"22 22 3 1",
". c None",
"# c #000000",
"a c #000082",
"......................",
".......#.....#........",
".......#.....#........",
".......#.....#........",
".......#....##........",
".......##...#.........",
"........#...#.........",
"........##.##.........",
".........###..........",
".........###..........",
"..........#...........",
".........a#a..........",
"........aa.aaa........",
".......a.a.a..a.......",
"......a..a.a...a......",
".....a...a.a....a.....",
"....a....a.a....a.....",
"....a....a..a...a.....",
"....a....a..a..a......",
"....a...a....aa.......",
".....aaa..............",
"......................"
]
image3_data = [
"22 22 6 1",
". c None",
"# c #000000",
"b c #000082",
"c c #3c3cfd",
"d c #8b8bfd",
"a c #ffffff",
"......................",
"......................",
"########..............",
"#aaaaaa##.............",
"#a####a#a#............",
"#aaaaaa#aa#...........",
"#a####a#bbbbbbbb......",
"#aaaaaa#baaaaaabb.....",
"#a#####aba####abcb....",
"#aaaaaaabaaaaaabdcb...",
"#a#####aba####abadcb..",
"#aaaaaaabaaaaaabbbbbb.",
"#a#####aba####aaaaaab.",
"#aaaaaaabaaaaaaaaaaab.",
"#a#####aba#########ab.",
"#aaaaaaabaaaaaaaaaaab.",
"########ba#########ab.",
"........baaaaaaaaaaab.",
"........ba#########ab.",
"........baaaaaaaaaaab.",
"........bbbbbbbbbbbbb.",
"......................"
]
image4_data = [
"22 22 8 1",
". c None",
"# c #000000",
"e c #000084",
"c c #848200",
"b c #848284",
"d c #c6c3c6",
"a c #ffff00",
"f c #ffffff",
"......................",
".......#####..........",
"..######aaa######.....",
".######aaaaa######....",
"##bcb##a###a##bcb##...",
"#bcb#ddddddddd#bcb#...",
"#cbc#ddddddddd#cbc#...",
"#bcb###########bcb#...",
"#cbcbcbcbcbcbcbcbc#...",
"#bcbcbcbcbcbcbcbcb#...",
"#cbcbcbceeeeeeeeee#...",
"#bcbcbcbefffffffefe...",
"#cbcbcbcefeeeeefeffe..",
"#bcbcbcbefffffffefffe.",
"#cbcbcbcefeeeeefeffffe",
"#bcbcbcbefffffffeeeeee",
"#cbcbcbcefeeeeeffffffe",
"#bcbcbcbeffffffffffffe",
"#cbcbcbcefeeeeeeeeeefe",
".#######effffffffffffe",
"........eeeeeeeeeeeeee",
"......................"
]

class MolDescriptorWin(QMainWindow):
  def __init__(self,parent = None,name = None,fl = 0):
    QMainWindow.__init__(self,parent,name,fl)
    self.statusBar()

    self.image0 = QPixmap(image0_data)
    self.image1 = QPixmap(image1_data)
    self.image2 = QPixmap(image2_data)
    self.image3 = QPixmap(image3_data)
    self.image4 = QPixmap(image4_data)

    if not name:
      self.setName("MolDescriptorWin")

    self.setIcon(self.image0)

    self.setCentralWidget(QWidget(self,"qt_central_widget"))
    MolDescriptorWinLayout = QVBoxLayout(self.centralWidget(),1,1,"MolDescriptorWinLayout")

    self.mainTab = QTabWidget(self.centralWidget(),"mainTab")

    self.sourcePage = QWidget(self.mainTab,"sourcePage")
    self.mainTab.insertTab(self.sourcePage,QString.fromLatin1(""))

    self.simplePage = QWidget(self.mainTab,"simplePage")
    self.mainTab.insertTab(self.simplePage,QString.fromLatin1(""))

    self.dataPage = QWidget(self.mainTab,"dataPage")
    self.mainTab.insertTab(self.dataPage,QString.fromLatin1(""))

    self.descPage = QWidget(self.mainTab,"descPage")
    self.mainTab.insertTab(self.descPage,QString.fromLatin1(""))
    MolDescriptorWinLayout.addWidget(self.mainTab)

    Layout8 = QHBoxLayout(None,0,6,"Layout8")
    Spacer7 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout8.addItem(Spacer7)

    self.buttonRead = QPushButton(self.centralWidget(),"buttonRead")
    self.buttonRead.setEnabled(0)
    Layout8.addWidget(self.buttonRead)
    Spacer6 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout8.addItem(Spacer6)

    self.buttonCalculate = QPushButton(self.centralWidget(),"buttonCalculate")
    self.buttonCalculate.setEnabled(0)
    Layout8.addWidget(self.buttonCalculate)
    Spacer8 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout8.addItem(Spacer8)
    MolDescriptorWinLayout.addLayout(Layout8)

    self.fileOpenAction = QAction(self,"fileOpenAction")
    self.fileOpenAction.setIconSet(QIconSet(self.image1))
    self.fileOpenAction.setAccel(4194383)
    self.fileSaveAsAction = QAction(self,"fileSaveAsAction")
    self.fileSaveAsAction.setAccel(0)
    self.fileCloseAction = QAction(self,"fileCloseAction")
    self.fileCloseAction.setAccel(0)
    self.editCutAction = QAction(self,"editCutAction")
    self.editCutAction.setIconSet(QIconSet(self.image2))
    self.editCutAction.setAccel(4194392)
    self.editCopyAction = QAction(self,"editCopyAction")
    self.editCopyAction.setIconSet(QIconSet(self.image3))
    self.editCopyAction.setAccel(4194371)
    self.editPasteAction = QAction(self,"editPasteAction")
    self.editPasteAction.setIconSet(QIconSet(self.image4))
    self.editPasteAction.setAccel(4194390)
    self.helpContentsAction = QAction(self,"helpContentsAction")
    self.helpContentsAction.setAccel(0)
    self.helpIndexAction = QAction(self,"helpIndexAction")
    self.helpIndexAction.setAccel(0)
    self.helpAboutAction = QAction(self,"helpAboutAction")
    self.helpAboutAction.setAccel(0)




    self.menubar = QMenuBar(self,"menubar")


    self.fileMenu = QPopupMenu(self)
    self.fileOpenAction.addTo(self.fileMenu)
    self.fileSaveAsAction.addTo(self.fileMenu)
    self.fileMenu.insertSeparator()
    self.fileMenu.insertSeparator()
    self.fileCloseAction.addTo(self.fileMenu)
    self.menubar.insertItem(QString(""),self.fileMenu,1)

    self.editMenu = QPopupMenu(self)
    self.editMenu.insertSeparator()
    self.editCutAction.addTo(self.editMenu)
    self.editCopyAction.addTo(self.editMenu)
    self.editPasteAction.addTo(self.editMenu)
    self.editMenu.insertSeparator()
    self.menubar.insertItem(QString(""),self.editMenu,2)

    self.helpMenu = QPopupMenu(self)
    self.helpContentsAction.addTo(self.helpMenu)
    self.helpIndexAction.addTo(self.helpMenu)
    self.helpMenu.insertSeparator()
    self.helpAboutAction.addTo(self.helpMenu)
    self.menubar.insertItem(QString(""),self.helpMenu,3)


    self.languageChange()

    self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.fileOpenAction,SIGNAL("activated()"),self.fileOpen)
    self.connect(self.fileSaveAsAction,SIGNAL("activated()"),self.fileSaveAs)
    self.connect(self.fileCloseAction,SIGNAL("activated()"),self.fileClose)
    self.connect(self.editCutAction,SIGNAL("activated()"),self.editCut)
    self.connect(self.editCopyAction,SIGNAL("activated()"),self.editCopy)
    self.connect(self.editPasteAction,SIGNAL("activated()"),self.editPaste)
    self.connect(self.helpIndexAction,SIGNAL("activated()"),self.helpIndex)
    self.connect(self.helpContentsAction,SIGNAL("activated()"),self.helpContents)
    self.connect(self.helpAboutAction,SIGNAL("activated()"),self.helpAbout)
    self.connect(self.buttonRead,SIGNAL("clicked()"),self.readClick)
    self.connect(self.buttonCalculate,SIGNAL("clicked()"),self.calcClick)


  def languageChange(self):
    self.setCaption(self.__tr("Molecular Descriptor Calculation"))
    self.mainTab.changeTab(self.sourcePage,self.__tr("Source"))
    self.mainTab.changeTab(self.simplePage,self.__tr("Simple Descriptors"))
    self.mainTab.changeTab(self.dataPage,self.__tr("Data"))
    self.mainTab.changeTab(self.descPage,self.__tr("Data+Descriptors"))
    self.buttonRead.setText(self.__tr("Read Compounds"))
    self.buttonCalculate.setText(self.__tr("Calculate Descriptors"))
    self.fileOpenAction.setText(self.__tr("Open"))
    self.fileOpenAction.setMenuText(self.__tr("&Open..."))
    self.fileSaveAsAction.setText(self.__tr("Save As"))
    self.fileSaveAsAction.setMenuText(self.__tr("Save &As..."))
    self.fileCloseAction.setText(self.__tr("Close"))
    self.fileCloseAction.setMenuText(self.__tr("C&lose"))
    self.editCutAction.setText(self.__tr("Cut"))
    self.editCutAction.setMenuText(self.__tr("&Cut"))
    self.editCopyAction.setText(self.__tr("Copy"))
    self.editCopyAction.setMenuText(self.__tr("C&opy"))
    self.editPasteAction.setText(self.__tr("Paste"))
    self.editPasteAction.setMenuText(self.__tr("&Paste"))
    self.helpContentsAction.setText(self.__tr("Contents"))
    self.helpContentsAction.setMenuText(self.__tr("&Contents..."))
    self.helpIndexAction.setText(self.__tr("Index"))
    self.helpIndexAction.setMenuText(self.__tr("&Index..."))
    self.helpAboutAction.setText(self.__tr("About"))
    self.helpAboutAction.setMenuText(self.__tr("&About..."))
    if self.menubar.findItem(1):
      self.menubar.findItem(1).setText(self.__tr("&File"))
    if self.menubar.findItem(2):
      self.menubar.findItem(2).setText(self.__tr("&Edit"))
    if self.menubar.findItem(3):
      self.menubar.findItem(3).setText(self.__tr("&Help"))


  def editCopy(self):
    print "MolDescriptorWin.editCopy(): Not implemented yet"

  def editCut(self):
    print "MolDescriptorWin.editCut(): Not implemented yet"

  def editFind(self):
    print "MolDescriptorWin.editFind(): Not implemented yet"

  def editPaste(self):
    print "MolDescriptorWin.editPaste(): Not implemented yet"

  def editRedo(self):
    print "MolDescriptorWin.editRedo(): Not implemented yet"

  def editUndo(self):
    print "MolDescriptorWin.editUndo(): Not implemented yet"

  def fileExit(self):
    print "MolDescriptorWin.fileExit(): Not implemented yet"

  def fileNew(self):
    print "MolDescriptorWin.fileNew(): Not implemented yet"

  def fileOpen(self):
    print "MolDescriptorWin.fileOpen(): Not implemented yet"

  def filePrint(self):
    print "MolDescriptorWin.filePrint(): Not implemented yet"

  def fileSave(self):
    print "MolDescriptorWin.fileSave(): Not implemented yet"

  def fileSaveAs(self):
    print "MolDescriptorWin.fileSaveAs(): Not implemented yet"

  def fileClose(self):
    print "MolDescriptorWin.fileClose(): Not implemented yet"

  def helpAbout(self):
    print "MolDescriptorWin.helpAbout(): Not implemented yet"

  def helpContents(self):
    print "MolDescriptorWin.helpContents(): Not implemented yet"

  def helpIndex(self):
    print "MolDescriptorWin.helpIndex(): Not implemented yet"

  def calcClick(self):
    print "MolDescriptorWin.calcClick(): Not implemented yet"

  def readClick(self):
    print "MolDescriptorWin.readClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("MolDescriptorWin",s,c)
