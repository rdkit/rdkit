# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SignatureViewDialog.ui'
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

class SignatureViewDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("SignatureViewDialog")

    self.setIcon(self.image0)

    SignatureViewDialogLayout = QVBoxLayout(self,2,2,"SignatureViewDialogLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")
    self.tabWidget.setEnabled(1)

    self.dbPage = QWidget(self.tabWidget,"dbPage")
    self.tabWidget.insertTab(self.dbPage,QString.fromLatin1(""))

    self.sigPage = QWidget(self.tabWidget,"sigPage")
    sigPageLayout = QVBoxLayout(self.sigPage,11,6,"sigPageLayout")

    Layout7 = QHBoxLayout(None,0,6,"Layout7")

    self.TextLabel1_2 = QLabel(self.sigPage,"TextLabel1_2")
    Layout7.addWidget(self.TextLabel1_2)

    self.sig_bitLine = QLineEdit(self.sigPage,"sig_bitLine")
    Layout7.addWidget(self.sig_bitLine)

    self.sig_buttonAdd = QPushButton(self.sigPage,"sig_buttonAdd")
    Layout7.addWidget(self.sig_buttonAdd)
    Spacer10 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout7.addItem(Spacer10)

    self.sig_buttonFile = QPushButton(self.sigPage,"sig_buttonFile")
    Layout7.addWidget(self.sig_buttonFile)
    sigPageLayout.addLayout(Layout7)
    self.tabWidget.insertTab(self.sigPage,QString.fromLatin1(""))
    SignatureViewDialogLayout.addWidget(self.tabWidget)

    Layout20 = QHBoxLayout(None,0,6,"Layout20")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout20.addWidget(self.buttonHelp)
    Spacer12 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer12)

    self.buttonRefresh = QPushButton(self,"buttonRefresh")
    self.buttonRefresh.setEnabled(0)
    Layout20.addWidget(self.buttonRefresh)
    Spacer12_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer12_2)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setAccel(0)
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout20.addWidget(self.buttonClose)
    SignatureViewDialogLayout.addLayout(Layout20)

    self.languageChange()

    self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.buttonClose,SIGNAL("clicked()"),self.accept)
    self.connect(self.buttonRefresh,SIGNAL("clicked()"),self.refreshClick)
    self.connect(self.sig_buttonAdd,SIGNAL("clicked()"),self.addBitClick)
    self.connect(self.sig_buttonFile,SIGNAL("clicked()"),self.loadBitsClick)


  def languageChange(self):
    self.setCaption(self.__tr("SignatureView"))
    self.tabWidget.changeTab(self.dbPage,self.__tr("Database"))
    self.TextLabel1_2.setText(self.__tr("Bit Number"))
    QToolTip.add(self.sig_bitLine,self.__tr("Enter the number of a bit to add to the worksheet."))
    self.sig_buttonAdd.setText(self.__tr("Add Bit"))
    self.sig_buttonFile.setText(self.__tr("Load Bits"))
    self.tabWidget.changeTab(self.sigPage,self.__tr("Signatures"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonRefresh.setText(self.__tr("Refresh"))
    QToolTip.add(self.buttonRefresh,self.__tr("Grabs signatures from the database"))
    self.buttonClose.setText(self.__tr("Close"))


  def addBitClick(self):
    print "SignatureViewDialog.addBitClick(): Not implemented yet"

  def loadBitsClick(self):
    print "SignatureViewDialog.loadBitsClick(): Not implemented yet"

  def refreshClick(self):
    print "SignatureViewDialog.refreshClick(): Not implemented yet"

  def screenClick(self):
    print "SignatureViewDialog.screenClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("SignatureViewDialog",s,c)
