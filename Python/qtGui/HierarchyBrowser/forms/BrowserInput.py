# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'forms/BrowserInput.ui'
#
# Created: Sat Sep 20 12:13:22 2003
#      by: The PyQt User Interface Compiler (pyuic) 3.7
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

class BrowserInput(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("BrowserInput")

    self.setIcon(self.image0)

    BrowserInputLayout = QGridLayout(self,1,1,2,4,"BrowserInputLayout")

    self.catalogLabel = QLabel(self,"catalogLabel")
    self.catalogLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    BrowserInputLayout.addWidget(self.catalogLabel,0,0)

    self.catalogEdit = QLineEdit(self,"catalogEdit")
    self.catalogEdit.setFrameShape(QLineEdit.LineEditPanel)
    self.catalogEdit.setFrameShadow(QLineEdit.Sunken)

    BrowserInputLayout.addMultiCellWidget(self.catalogEdit,0,0,1,2)

    self.catalogBrowseButton = QToolButton(self,"catalogBrowseButton")

    BrowserInputLayout.addWidget(self.catalogBrowseButton,0,3)

    self.gainsBrowseButton = QToolButton(self,"gainsBrowseButton")

    BrowserInputLayout.addWidget(self.gainsBrowseButton,1,3)

    self.gainsEdit = QLineEdit(self,"gainsEdit")

    BrowserInputLayout.addMultiCellWidget(self.gainsEdit,1,1,1,2)

    self.catalogLabel_2 = QLabel(self,"catalogLabel_2")
    self.catalogLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    BrowserInputLayout.addWidget(self.catalogLabel_2,1,0)
    spacer = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    BrowserInputLayout.addItem(spacer,2,2)

    self.catalogLabel_2_2 = QLabel(self,"catalogLabel_2_2")
    self.catalogLabel_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    BrowserInputLayout.addWidget(self.catalogLabel_2_2,2,0)

    self.catalogLabel_2_2_2 = QLabel(self,"catalogLabel_2_2_2")
    self.catalogLabel_2_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    BrowserInputLayout.addWidget(self.catalogLabel_2_2_2,3,0)
    spacer_2 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    BrowserInputLayout.addItem(spacer_2,3,2)
    spacer_3 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    BrowserInputLayout.addItem(spacer_3,4,2)

    self.catalogLabel_2_2_2_2 = QLabel(self,"catalogLabel_2_2_2_2")
    self.catalogLabel_2_2_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    BrowserInputLayout.addWidget(self.catalogLabel_2_2_2_2,4,0)

    self.maxLevelSpin = QSpinBox(self,"maxLevelSpin")
    self.maxLevelSpin.setEnabled(0)

    BrowserInputLayout.addWidget(self.maxLevelSpin,4,1)

    self.minLevelSpin = QSpinBox(self,"minLevelSpin")
    self.minLevelSpin.setEnabled(0)

    BrowserInputLayout.addWidget(self.minLevelSpin,3,1)

    self.numBitsSpin = QSpinBox(self,"numBitsSpin")
    self.numBitsSpin.setEnabled(0)
    self.numBitsSpin.setMaxValue(200)
    self.numBitsSpin.setLineStep(10)
    self.numBitsSpin.setValue(50)

    BrowserInputLayout.addWidget(self.numBitsSpin,2,1)

    self.languageChange()

    self.resize(QSize(296,126).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.catalogBrowseButton,SIGNAL("clicked()"),self.catalogBrowseClicked)
    self.connect(self.gainsBrowseButton,SIGNAL("clicked()"),self.gainsBrowseClicked)
    self.connect(self.catalogEdit,SIGNAL("returnPressed()"),self.catalogReturnPressed)
    self.connect(self.gainsEdit,SIGNAL("returnPressed()"),self.gainsReturnPressed)


  def languageChange(self):
    self.setCaption(self.__tr("HierarchyBrowserInput"))
    self.catalogLabel.setText(self.__tr("Catalog"))
    QToolTip.add(self.catalogEdit,self.__tr("Enter the name of the catalog file here."))
    self.catalogBrowseButton.setText(self.__tr("..."))
    self.gainsBrowseButton.setText(self.__tr("..."))
    QToolTip.add(self.gainsEdit,self.__tr("Enter the name of the gains file here."))
    self.catalogLabel_2.setText(self.__tr("Gains"))
    self.catalogLabel_2_2.setText(self.__tr("Num Bits"))
    QToolTip.add(self.catalogLabel_2_2,QString.null)
    self.catalogLabel_2_2_2.setText(self.__tr("Min Level"))
    QToolTip.add(self.catalogLabel_2_2_2,QString.null)
    self.catalogLabel_2_2_2_2.setText(self.__tr("Max Level"))
    QToolTip.add(self.catalogLabel_2_2_2_2,QString.null)
    QToolTip.add(self.maxLevelSpin,self.__tr("The maximum hierarchy level (e.g. path length) to be included in the browser."))
    QToolTip.add(self.minLevelSpin,self.__tr("The minimum hierarchy level (e.g. path length) to be included in the browser."))
    QToolTip.add(self.numBitsSpin,self.__tr("The number of top bits to be included in the browser."))


  def catalogBrowseClicked(self):
    print "BrowserInput.catalogBrowseClicked(): Not implemented yet"

  def gainsBrowseClicked(self):
    print "BrowserInput.gainsBrowseClicked(): Not implemented yet"

  def catalogReturnPressed(self):
    print "BrowserInput.catalogReturnPressed(): Not implemented yet"

  def gainsReturnPressed(self):
    print "BrowserInput.gainsReturnPressed(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("BrowserInput",s,c)
