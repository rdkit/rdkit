# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/CompositeWindow.ui'
#
# Created: Sun Apr 23 08:07:30 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.11
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
"22 22 7 1",
". c None",
"# c #000000",
"b c #2e2e2e",
"c c #5c5c5c",
"d c #878787",
"e c #c2c2c2",
"a c #ffffff",
"......................",
"....##########........",
"....#aaaaaaa#b#.......",
"....#aaaaaaa#cb#......",
"....#aaaaaaa#dcb#.....",
"....#aaaaaaa#edcb#....",
"....#aaaaaaa#aedcb#...",
"....#aaaaaaa#######...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....#aaaaaaaaaaaaa#...",
"....###############...",
"......................",
"......................"
]
image2_data = [
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
image3_data = [
"22 22 88 2",
"Qt c None",
".2 c #000000",
".S c #08ff08",
"#v c #100810",
".U c #101010",
"#c c #101018",
".M c #181018",
"#e c #181818",
".A c #181821",
".L c #211821",
"#l c #212121",
".z c #212129",
"#m c #292129",
"#u c #292929",
"#n c #292931",
".R c #29ff29",
"#o c #312931",
".T c #313131",
"#p c #313139",
".Z c #31ff31",
"#q c #393139",
"#t c #393939",
".y c #393942",
"#s c #423942",
".o c #424242",
"#h c #4a4a52",
".n c #5a525a",
"#r c #5a5a63",
".I c #5ace5a",
"#b c #6b636b",
".p c #6b6b6b",
".x c #6b6b73",
".Y c #6bff63",
".l c #736b73",
".t c #7b737b",
".s c #7b7384",
".0 c #7bff7b",
".r c #847b84",
".u c #847b8c",
"#g c #84848c",
".v c #8c7b94",
"#i c #8c848c",
".w c #8c8494",
"#j c #8c8c8c",
".8 c #8c8c94",
".m c #948c94",
"#k c #948c9c",
"#f c #949494",
".q c #94949c",
".J c #94c694",
"#d c #9c949c",
"#a c #9c94a5",
".k c #9c9c9c",
".N c #9c9ca5",
".H c #9ccea5",
".K c #a59ca5",
"#. c #a59cad",
".i c #a5a5a5",
".3 c #a5a5ad",
"## c #ad9cad",
".V c #ada5ad",
".d c #adadad",
".j c #adadb5",
".9 c #b5adb5",
".# c #b5b5b5",
".a c #bdbdbd",
".7 c #bdd6bd",
".c c #c6c6c6",
".5 c #cec6ce",
".b c #cecece",
".4 c #ceced6",
".F c #d6ced6",
".G c #d6cede",
".h c #d6d6d6",
".E c #d6d6de",
".Q c #d6ffd6",
".B c #ded6de",
".1 c #ded6e7",
".g c #dedede",
".D c #dedee7",
".6 c #e7dee7",
".f c #e7e7e7",
".C c #e7e7ef",
".X c #e7ffe7",
".O c #efe7ef",
".e c #efefef",
".W c #f7f7f7",
".P c #ffffff",
"QtQtQtQtQtQt.#.a.b.b.b.b.c.c.a.a.d.aQtQtQtQt",
"QtQtQtQtQtQt.a.e.f.f.f.f.f.e.e.e.g.aQtQtQtQt",
"QtQtQtQtQtQt.a.c.c.c.b.b.c.c.c.c.a.cQtQtQtQt",
"QtQtQtQtQtQt.#.a.a.a.a.#.a.a.#.#.d.aQtQtQtQt",
"QtQtQtQtQt.c.d.c.a.c.c.c.a.a.a.c.#QtQtQtQtQt",
"QtQtQtQtQt.a.a.#.a.a.a.a.a.a.c.c.#QtQtQtQtQt",
"QtQtQtQtQt.a.#.c.a.a.a.a.a.c.a.c.dQtQtQtQtQt",
"QtQtQtQtQt.c.a.a.a.a.a.a.a.a.a.a.#QtQtQtQtQt",
"QtQtQtQtQt.d.b.f.g.g.g.g.g.g.h.g.i.i.jQtQtQt",
"QtQtQt.a.k.l.#.h.b.h.b.h.b.h.g.g.m.n.o.p.#Qt",
"QtQt.a.q.r.s.t.t.t.t.t.t.t.u.v.w.x.y.z.A.o.i",
"Qt.a.k.B.C.D.B.E.E.E.E.F.G.H.I.J.K.o.L.L.M.y",
".a.N.O.P.P.P.P.P.P.P.P.P.Q.R.S.R.b.v.T.A.U.L",
".V.W.P.P.P.P.P.P.P.P.P.P.X.Y.Z.0.P.1.t.A.2.L",
".3.E.4.5.4.h.E.E.g.6.D.B.D.E.7.F.4.5.8.M.2.A",
".m.9.j.V.3#..3.K#.#..i#..K#.###a.q.8#b#c.2.L",
".m.j.j#..3.K.K.K.N.K.N.N.N.N#a#d#d.w#b#c.2#e",
"#f#.#..K.N.K.N.N.N#a.k#a#d#d#d#a.m#g#b.M.2#h",
".m.3.K.K#a.k#a#d#a.k#a#d#a#d.q.m.8#i.x#c#e.d",
"#f#g#i.w#j.w#i.8.w#i.8.8.m.8.m#k.8.w#b#e#fQt",
".#.l.z.A#l.z#m#m#m#n#o#o#p#p#q#q#p#o#p#fQtQt",
"QtQt.d#r#s#s#t#p.T.T.T#u#u.z#e#e#v.o.kQtQtQt"
]

class CompositeWindow(QMainWindow):
  def __init__(self,parent = None,name = None,fl = 0):
    QMainWindow.__init__(self,parent,name,fl)
    self.statusBar()

    self.image0 = QPixmap(image0_data)
    self.image1 = QPixmap(image1_data)
    self.image2 = QPixmap(image2_data)
    self.image3 = QPixmap(image3_data)

    if not name:
      self.setName("CompositeWindow")

    self.setIcon(self.image0)

    self.setCentralWidget(QWidget(self,"qt_central_widget"))
    CompositeWindowLayout = QVBoxLayout(self.centralWidget(),2,2,"CompositeWindowLayout")

    self.tabWidget = QTabWidget(self.centralWidget(),"tabWidget")
    CompositeWindowLayout.addWidget(self.tabWidget)

    self.fileNewAction = QAction(self,"fileNewAction")
    self.fileNewAction.setIconSet(QIconSet(self.image1))
    self.fileNewAction.setAccel(4194382)
    self.fileOpenAction = QAction(self,"fileOpenAction")
    self.fileOpenAction.setIconSet(QIconSet(self.image2))
    self.fileOpenAction.setAccel(4194383)
    self.fileSaveAsAction = QAction(self,"fileSaveAsAction")
    self.fileSaveAsAction.setAccel(0)
    self.filePrintAction = QAction(self,"filePrintAction")
    self.filePrintAction.setIconSet(QIconSet(self.image3))
    self.filePrintAction.setAccel(4194384)
    self.fileCloseAction = QAction(self,"fileCloseAction")
    self.fileCloseAction.setAccel(0)
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
    self.filePrintAction.addTo(self.fileMenu)
    self.fileMenu.insertSeparator()
    self.fileCloseAction.addTo(self.fileMenu)
    self.menubar.insertItem(QString(""),self.fileMenu,1)

    self.helpMenu = QPopupMenu(self)
    self.helpContentsAction.addTo(self.helpMenu)
    self.helpMenu.insertSeparator()
    self.helpAboutAction.addTo(self.helpMenu)
    self.menubar.insertItem(QString(""),self.helpMenu,2)


    self.languageChange()

    self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.fileNewAction,SIGNAL("activated()"),self.fileNew)
    self.connect(self.fileOpenAction,SIGNAL("activated()"),self.fileOpen)
    self.connect(self.fileSaveAsAction,SIGNAL("activated()"),self.fileSaveAs)
    self.connect(self.filePrintAction,SIGNAL("activated()"),self.filePrint)
    self.connect(self.fileCloseAction,SIGNAL("activated()"),self.fileClose)
    self.connect(self.helpIndexAction,SIGNAL("activated()"),self.helpIndex)
    self.connect(self.helpContentsAction,SIGNAL("activated()"),self.helpContents)
    self.connect(self.helpAboutAction,SIGNAL("activated()"),self.helpAbout)


  def languageChange(self):
    self.setCaption(self.__tr("Composite"))
    self.fileNewAction.setText(self.__tr("New"))
    self.fileNewAction.setMenuText(self.__tr("&New"))
    self.fileOpenAction.setText(self.__tr("Open"))
    self.fileOpenAction.setMenuText(self.__tr("&Open..."))
    self.fileSaveAsAction.setText(self.__tr("Save As"))
    self.fileSaveAsAction.setMenuText(self.__tr("Save &As..."))
    self.filePrintAction.setText(self.__tr("Print"))
    self.filePrintAction.setMenuText(self.__tr("&Print..."))
    self.fileCloseAction.setText(self.__tr("Close"))
    self.fileCloseAction.setMenuText(self.__tr("&Close"))
    self.helpContentsAction.setText(self.__tr("Contents"))
    self.helpContentsAction.setMenuText(self.__tr("&Contents..."))
    self.helpIndexAction.setText(self.__tr("Index"))
    self.helpIndexAction.setMenuText(self.__tr("&Index..."))
    self.helpAboutAction.setText(self.__tr("About"))
    self.helpAboutAction.setMenuText(self.__tr("&About..."))
    if self.menubar.findItem(1):
      self.menubar.findItem(1).setText(self.__tr("&File"))
    if self.menubar.findItem(2):
      self.menubar.findItem(2).setText(self.__tr("&Help"))


  def fileNew(self):
    print "CompositeWindow.fileNew(): Not implemented yet"

  def fileOpen(self):
    print "CompositeWindow.fileOpen(): Not implemented yet"

  def fileSaveAs(self):
    print "CompositeWindow.fileSaveAs(): Not implemented yet"

  def filePrint(self):
    print "CompositeWindow.filePrint(): Not implemented yet"

  def fileClose(self):
    print "CompositeWindow.fileClose(): Not implemented yet"

  def helpIndex(self):
    print "CompositeWindow.helpIndex(): Not implemented yet"

  def helpContents(self):
    print "CompositeWindow.helpContents(): Not implemented yet"

  def helpAbout(self):
    print "CompositeWindow.helpAbout(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("CompositeWindow",s,c)
