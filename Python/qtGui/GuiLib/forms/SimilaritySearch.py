# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SimilaritySearch.ui'
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

class SimilaritySearchWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("SimilaritySearchWidget")

    self.setIcon(self.image0)

    SimilaritySearchWidgetLayout = QVBoxLayout(self,4,2,"SimilaritySearchWidgetLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")

    self.queryPage = QWidget(self.tabWidget,"queryPage")
    self.tabWidget.insertTab(self.queryPage,QString.fromLatin1(""))

    self.databasePage = QWidget(self.tabWidget,"databasePage")
    self.tabWidget.insertTab(self.databasePage,QString.fromLatin1(""))

    self.paramsPage = QWidget(self.tabWidget,"paramsPage")
    self.tabWidget.insertTab(self.paramsPage,QString.fromLatin1(""))
    SimilaritySearchWidgetLayout.addWidget(self.tabWidget)

    layout2 = QHBoxLayout(None,0,2,"layout2")

    self.searchButton = QPushButton(self,"searchButton")
    self.searchButton.setEnabled(0)
    layout2.addWidget(self.searchButton)
    spacer1 = QSpacerItem(338,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout2.addItem(spacer1)
    SimilaritySearchWidgetLayout.addLayout(layout2)

    self.languageChange()

    self.resize(QSize(444,219).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.searchButton,SIGNAL("clicked()"),self.searchButtonClicked)


  def languageChange(self):
    self.setCaption(self.__tr("SimilaritySearch"))
    self.tabWidget.changeTab(self.queryPage,self.__tr("&Query"))
    self.tabWidget.changeTab(self.databasePage,self.__tr("&Database"))
    self.tabWidget.changeTab(self.paramsPage,self.__tr("&Parameters"))
    self.searchButton.setText(self.__tr("Search"))
    QToolTip.add(self.searchButton,self.__tr("Carry out the similarity search."))


  def searchButtonClicked(self):
    print "SimilaritySearchWidget.searchButtonClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("SimilaritySearchWidget",s,c)
