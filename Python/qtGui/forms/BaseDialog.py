# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'BaseDialog.ui'
#
# Created: Mon Jan 23 08:34:47 2006
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

class BaseDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("BaseDialog")

    self.setIcon(self.image0)
    self.setSizeGripEnabled(0)


    LayoutWidget = QWidget(self,"Layout1")
    LayoutWidget.setGeometry(QRect(11,310,596,32))
    Layout1 = QHBoxLayout(LayoutWidget,0,6,"Layout1")

    self.buttonHelp = QPushButton(LayoutWidget,"buttonHelp")
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    Horizontal_Spacing2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(Horizontal_Spacing2)

    self.buttonOk = QPushButton(LayoutWidget,"buttonOk")
    self.buttonOk.setAccel(0)
    self.buttonOk.setAutoDefault(1)
    self.buttonOk.setDefault(1)
    Layout1.addWidget(self.buttonOk)

    self.buttonCancel = QPushButton(LayoutWidget,"buttonCancel")
    self.buttonCancel.setAccel(0)
    self.buttonCancel.setAutoDefault(1)
    Layout1.addWidget(self.buttonCancel)

    self.languageChange()

    self.resize(QSize(618,352).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self.accept)
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self.reject)


  def languageChange(self):
    self.setCaption(self.__tr("Dialog"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonOk.setText(self.__tr("OK"))
    self.buttonCancel.setText(self.__tr("Cancel"))


  def dbFileClick(self):
    print "BaseDialog.dbFileClick(): Not implemented yet"

  def fileClick(self):
    print "BaseDialog.fileClick(): Not implemented yet"

  def refreshContents(self):
    print "BaseDialog.refreshContents(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("BaseDialog",s,c)
