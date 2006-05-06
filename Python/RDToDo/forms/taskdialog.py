# Form implementation generated from reading ui file 'taskdialog.ui'
#
# Created: Fri Feb 21 07:26:53 2003
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class TaskEdit(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    if name == None:
      self.setName("taskEdit")

    self.resize(542,335)
    self.setCaption(self.trUtf8("Edit Task"))
    self.setSizeGripEnabled(1)

    taskEditLayout = QVBoxLayout(self,4,4,"taskEditLayout")

    self.notebook = QTabWidget(self,"notebook")

    self.tab = QWidget(self.notebook,"tab")
    tabLayout = QVBoxLayout(self.tab,4,4,"tabLayout")

    Layout8 = QHBoxLayout(None,0,6,"Layout8")

    self.TextLabel1 = QLabel(self.tab,"TextLabel1")
    self.TextLabel1.setText(self.trUtf8("Name"))
    Layout8.addWidget(self.TextLabel1)

    self.nameEdit = QLineEdit(self.tab,"nameEdit")
    Layout8.addWidget(self.nameEdit)
    tabLayout.addLayout(Layout8)

    Layout51 = QHBoxLayout(None,0,6,"Layout51")

    self.dueButton = QCheckBox(self.tab,"dueButton")
    self.dueButton.setText(self.trUtf8("Due?"))
    QToolTip.add(self.dueButton,self.trUtf8("toggles use of the due date"))
    Layout51.addWidget(self.dueButton)

    self.dueBox = QDateEdit(self.tab,"dueBox")
    self.dueBox.setEnabled(0)
    QToolTip.add(self.dueBox,self.trUtf8("the date the task is due."))
    Layout51.addWidget(self.dueBox)

    self.priorityButton = QCheckBox(self.tab,"priorityButton")
    self.priorityButton.setText(self.trUtf8("Priority?"))
    QToolTip.add(self.priorityButton,self.trUtf8("priority of the item on a scale of 0-99"))
    Layout51.addWidget(self.priorityButton)

    self.priorityBox = QSpinBox(self.tab,"priorityBox")
    self.priorityBox.setEnabled(0)
    QToolTip.add(self.priorityBox,self.trUtf8("priority of the item on a scale of 0-99"))
    Layout51.addWidget(self.priorityBox)
    spacer = QSpacerItem(60,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout51.addItem(spacer)
    tabLayout.addLayout(Layout51)

    self.TextLabel4 = QLabel(self.tab,"TextLabel4")
    self.TextLabel4.setText(self.trUtf8("Note:"))
    tabLayout.addWidget(self.TextLabel4)

    self.note_edit = QTextEdit(self.tab,"note_edit")
    tabLayout.addWidget(self.note_edit)
    self.notebook.insertTab(self.tab,self.trUtf8("Info"))
    taskEditLayout.addWidget(self.notebook)

    Layout1 = QHBoxLayout(None,0,6,"Layout1")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setEnabled(0)
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer_2)

    self.buttonOk = QPushButton(self,"buttonOk")
    self.buttonOk.setText(self.trUtf8("OK"))
    self.buttonOk.setAccel(0)
    self.buttonOk.setAutoDefault(1)
    self.buttonOk.setDefault(1)
    Layout1.addWidget(self.buttonOk)

    self.buttonCancel = QPushButton(self,"buttonCancel")
    self.buttonCancel.setText(self.trUtf8("Cancel"))
    self.buttonCancel.setAccel(0)
    self.buttonCancel.setAutoDefault(1)
    Layout1.addWidget(self.buttonCancel)
    taskEditLayout.addLayout(Layout1)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.priorityButton,SIGNAL("clicked()"),self.priorityClick)
    self.connect(self.dueButton,SIGNAL("clicked()"),self.dueClick)

    self.setTabOrder(self.nameEdit,self.dueButton)
    self.setTabOrder(self.dueButton,self.dueBox)
    self.setTabOrder(self.dueBox,self.priorityButton)
    self.setTabOrder(self.priorityButton,self.priorityBox)
    self.setTabOrder(self.priorityBox,self.note_edit)
    self.setTabOrder(self.note_edit,self.buttonOk)
    self.setTabOrder(self.buttonOk,self.buttonCancel)
    self.setTabOrder(self.buttonCancel,self.buttonHelp)
    self.setTabOrder(self.buttonHelp,self.notebook)

  def dueClick(self):
    print "TaskEdit.dueClick(): Not implemented yet"

  def priorityClick(self):
    print "TaskEdit.priorityClick(): Not implemented yet"
