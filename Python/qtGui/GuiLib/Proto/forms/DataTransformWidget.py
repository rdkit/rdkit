# Form implementation generated from reading ui file 'forms/DataTransformWidget.ui'
#
# Created: Wed Nov 27 21:50:51 2002
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *
from qttable import QTable

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

class DataTransformWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    image0 = QPixmap(image0_data)

    if name == None:
      self.setName("DataTransformWidget")

    self.resize(600,162)
    self.setCaption(self.trUtf8("Data Transformation"))
    self.setIcon(image0)

    DataTransformWidgetLayout = QVBoxLayout(self,2,4,"DataTransformWidgetLayout")

    self.choiceTable = QTable(self,"choiceTable")
    self.choiceTable.setNumCols(self.choiceTable.numCols() + 1)
    self.choiceTable.horizontalHeader().setLabel(self.choiceTable.numCols() - 1,self.trUtf8("Name"))
    self.choiceTable.setNumCols(self.choiceTable.numCols() + 1)
    self.choiceTable.horizontalHeader().setLabel(self.choiceTable.numCols() - 1,self.trUtf8("Apply To"))
    self.choiceTable.setNumCols(self.choiceTable.numCols() + 1)
    self.choiceTable.horizontalHeader().setLabel(self.choiceTable.numCols() - 1,self.trUtf8("Operation"))
    self.choiceTable.setNumCols(self.choiceTable.numCols() + 1)
    self.choiceTable.horizontalHeader().setLabel(self.choiceTable.numCols() - 1,self.trUtf8("Formula"))
    self.choiceTable.setNumRows(self.choiceTable.numRows() + 1)
    self.choiceTable.verticalHeader().setLabel(self.choiceTable.numRows() - 1,self.trUtf8("0"))
    self.choiceTable.setNumRows(self.choiceTable.numRows() + 1)
    self.choiceTable.verticalHeader().setLabel(self.choiceTable.numRows() - 1,self.trUtf8("1"))
    self.choiceTable.setNumRows(self.choiceTable.numRows() + 1)
    self.choiceTable.verticalHeader().setLabel(self.choiceTable.numRows() - 1,self.trUtf8("2"))
    self.choiceTable.setResizePolicy(QTable.Default)
    self.choiceTable.setNumRows(3)
    self.choiceTable.setNumCols(4)
    self.choiceTable.setRowMovingEnabled(1)
    self.choiceTable.setSelectionMode(QTable.NoSelection)
    self.choiceTable.setFocusStyle(QTable.FollowStyle)
    DataTransformWidgetLayout.addWidget(self.choiceTable)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    Layout1.addWidget(self.buttonHelp)

    self.buttonAdd = QPushButton(self,"buttonAdd")
    self.buttonAdd.setText(self.trUtf8("Add New Operation"))
    QToolTip.add(self.buttonAdd,self.trUtf8("Add a new row (operation) to the transform"))
    Layout1.addWidget(self.buttonAdd)
    spacer = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer)

    self.buttonSave = QPushButton(self,"buttonSave")
    self.buttonSave.setText(self.trUtf8("Save"))
    QToolTip.add(self.buttonSave,self.trUtf8("Save this data transformation object"))
    Layout1.addWidget(self.buttonSave)

    self.buttonLoad = QPushButton(self,"buttonLoad")
    self.buttonLoad.setText(self.trUtf8("Load"))
    QToolTip.add(self.buttonLoad,self.trUtf8("Load a data transformation object"))
    Layout1.addWidget(self.buttonLoad)
    DataTransformWidgetLayout.addLayout(Layout1)

    self.connect(self.buttonHelp,SIGNAL("clicked()"),self.helpClick)
    self.connect(self.buttonAdd,SIGNAL("clicked()"),self.addClick)
    self.connect(self.buttonSave,SIGNAL("pressed()"),self.saveClick)
    self.connect(self.buttonLoad,SIGNAL("clicked()"),self.loadClick)

  def addClick(self):
    print "DataTransformWidget.addClick(): Not implemented yet"

  def helpClick(self):
    print "DataTransformWidget.helpClick(): Not implemented yet"

  def loadClick(self):
    print "DataTransformWidget.loadClick(): Not implemented yet"

  def saveClick(self):
    print "DataTransformWidget.saveClick(): Not implemented yet"
