# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'PubmedSearchWidget.ui'
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

class PubmedSearchWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("PubmedSearchWidget")

    self.setIcon(self.image0)

    PubmedSearchWidgetLayout = QVBoxLayout(self,4,2,"PubmedSearchWidgetLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")

    self.queryPage = QWidget(self.tabWidget,"queryPage")
    queryPageLayout = QVBoxLayout(self.queryPage,4,2,"queryPageLayout")

    layout2 = QHBoxLayout(None,0,2,"layout2")

    self.textLabel2 = QLabel(self.queryPage,"textLabel2")
    self.textLabel2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    layout2.addWidget(self.textLabel2)

    self.queryLine = QLineEdit(self.queryPage,"queryLine")
    layout2.addWidget(self.queryLine)
    queryPageLayout.addLayout(layout2)
    self.tabWidget.insertTab(self.queryPage,QString.fromLatin1(""))
    PubmedSearchWidgetLayout.addWidget(self.tabWidget)

    layout3 = QHBoxLayout(None,0,2,"layout3")

    self.searchButton = QPushButton(self,"searchButton")
    self.searchButton.setEnabled(1)
    layout3.addWidget(self.searchButton)
    spacer1 = QSpacerItem(30,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout3.addItem(spacer1)

    self.recordsButton = QPushButton(self,"recordsButton")
    self.recordsButton.setEnabled(0)
    layout3.addWidget(self.recordsButton)
    spacer2 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout3.addItem(spacer2)

    self.relatedButton = QPushButton(self,"relatedButton")
    self.relatedButton.setEnabled(0)
    layout3.addWidget(self.relatedButton)
    spacer2_2 = QSpacerItem(30,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout3.addItem(spacer2_2)

    self.exportButton = QPushButton(self,"exportButton")
    self.exportButton.setEnabled(0)
    layout3.addWidget(self.exportButton)
    PubmedSearchWidgetLayout.addLayout(layout3)

    self.languageChange()

    self.resize(QSize(491,345).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.searchButton,SIGNAL("clicked()"),self.searchButtonClicked)
    self.connect(self.recordsButton,SIGNAL("clicked()"),self.recordsButtonClicked)
    self.connect(self.relatedButton,SIGNAL("clicked()"),self.relatedButtonClicked)
    self.connect(self.exportButton,SIGNAL("clicked()"),self.exportButtonClicked)


  def languageChange(self):
    self.setCaption(self.__tr("PubmedSearch"))
    self.textLabel2.setText(self.__tr("Query:"))
    QToolTip.add(self.textLabel2,self.__tr("This is the actual query that will be sent to Pubmed"))
    QToolTip.add(self.queryLine,self.__tr("This is the actual query that will be sent to Pubmed"))
    self.tabWidget.changeTab(self.queryPage,self.__tr("Query"))
    self.searchButton.setText(self.__tr("Search"))
    QToolTip.add(self.searchButton,self.__tr("Search Pubmed for the current query."))
    self.recordsButton.setText(self.__tr("Records"))
    QToolTip.add(self.recordsButton,self.__tr("Retrieve full records for the active result set."))
    self.relatedButton.setText(self.__tr("Related"))
    QToolTip.add(self.relatedButton,self.__tr("Retrieve information about related records."))
    self.exportButton.setText(self.__tr("Export"))
    QToolTip.add(self.exportButton,self.__tr("Exports the selected records to a text file readable by EndNote."))


  def searchButtonClicked(self):
    print "PubmedSearchWidget.searchButtonClicked(): Not implemented yet"

  def recordsButtonClicked(self):
    print "PubmedSearchWidget.recordsButtonClicked(): Not implemented yet"

  def relatedButtonClicked(self):
    print "PubmedSearchWidget.relatedButtonClicked(): Not implemented yet"

  def exportButtonClicked(self):
    print "PubmedSearchWidget.exportButtonClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("PubmedSearchWidget",s,c)
