# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ModelBrowserDialog.ui'
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

class ModelBrowserDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("ModelBrowserDialog")

    self.setIcon(self.image0)

    ModelBrowserDialogLayout = QVBoxLayout(self,2,2,"ModelBrowserDialogLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")
    self.tabWidget.setEnabled(1)

    self.dbPage = QWidget(self.tabWidget,"dbPage")
    self.tabWidget.insertTab(self.dbPage,QString.fromLatin1(""))

    self.modelsPage = QWidget(self.tabWidget,"modelsPage")
    modelsPageLayout = QVBoxLayout(self.modelsPage,2,2,"modelsPageLayout")
    self.tabWidget.insertTab(self.modelsPage,QString.fromLatin1(""))
    ModelBrowserDialogLayout.addWidget(self.tabWidget)

    Layout7 = QHBoxLayout(None,0,6,"Layout7")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setAutoDefault(1)
    Layout7.addWidget(self.buttonHelp)
    Spacer12 = QSpacerItem(62,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout7.addItem(Spacer12)

    self.buttonRefresh = QPushButton(self,"buttonRefresh")
    self.buttonRefresh.setEnabled(0)
    Layout7.addWidget(self.buttonRefresh)
    Spacer13 = QSpacerItem(62,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout7.addItem(Spacer13)

    self.buttonTransfer = QPushButton(self,"buttonTransfer")
    self.buttonTransfer.setEnabled(0)
    Layout7.addWidget(self.buttonTransfer)
    Spacer15_2 = QSpacerItem(90,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout7.addItem(Spacer15_2)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout7.addWidget(self.buttonClose)
    ModelBrowserDialogLayout.addLayout(Layout7)

    self.languageChange()

    self.resize(QSize(430,262).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.buttonClose,SIGNAL("clicked()"),self.accept)
    self.connect(self.buttonTransfer,SIGNAL("clicked()"),self.transferClick)
    self.connect(self.buttonRefresh,SIGNAL("clicked()"),self.refreshClick)


  def languageChange(self):
    self.setCaption(self.__tr("ModelBrowser"))
    self.tabWidget.changeTab(self.dbPage,self.__tr("Database"))
    self.tabWidget.changeTab(self.modelsPage,self.__tr("Models"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonHelp.setAccel(self.__tr("F1"))
    self.buttonRefresh.setText(self.__tr("Refresh"))
    QToolTip.add(self.buttonRefresh,self.__tr("Updates the contents of the Models tab from the database."))
    self.buttonTransfer.setText(self.__tr("Transfer"))
    QToolTip.add(self.buttonTransfer,self.__tr("Transfer the selected models to the Composite Interaction window"))
    self.buttonClose.setText(self.__tr("Close"))
    self.buttonClose.setAccel(QString.null)


  def transferClick(self):
    print "ModelBrowserDialog.transferClick(): Not implemented yet"

  def refreshClick(self):
    print "ModelBrowserDialog.refreshClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("ModelBrowserDialog",s,c)
