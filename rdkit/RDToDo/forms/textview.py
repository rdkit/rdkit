# Form implementation generated from reading ui file 'forms/textview.ui'
#
# Created: Mon Apr 15 16:55:46 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class TextView(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if name == None:
      self.setName("TextView")

    self.resize(600,480)
    self.setCaption(self.trUtf8("Text View"))

    TextViewLayout = QHBoxLayout(self,11,6,"TextViewLayout")

    self.textBrowser = QTextEdit(self,"textBrowser")
    TextViewLayout.addWidget(self.textBrowser)
