# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ScreenCompositeDialog.ui'
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

class ScreenCompositeDialog(QDialog):
  def __init__(self,parent = None,name = None,modal = 0,fl = 0):
    QDialog.__init__(self,parent,name,modal,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("ScreenCompositeDialog")

    self.setIcon(self.image0)

    ScreenCompositeDialogLayout = QVBoxLayout(self,4,4,"ScreenCompositeDialogLayout")

    self.tabWidget = QTabWidget(self,"tabWidget")
    self.tabWidget.setEnabled(1)

    self.optionsPage = QWidget(self.tabWidget,"optionsPage")
    optionsPageLayout = QVBoxLayout(self.optionsPage,2,2,"optionsPageLayout")

    Layout18 = QHBoxLayout(None,0,6,"Layout18")

    self.options_thresholdEdit = QLineEdit(self.optionsPage,"options_thresholdEdit")
    self.options_thresholdEdit.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.options_thresholdEdit.sizePolicy().hasHeightForWidth()))
    self.options_thresholdEdit.setMinimumSize(QSize(50,0))
    self.options_thresholdEdit.setMaximumSize(QSize(50,32767))
    Layout18.addWidget(self.options_thresholdEdit)

    self.TextLabel1 = QLabel(self.optionsPage,"TextLabel1")
    Layout18.addWidget(self.TextLabel1)
    Spacer21 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout18.addItem(Spacer21)
    optionsPageLayout.addLayout(Layout18)

    Layout19 = QHBoxLayout(None,0,6,"Layout19")

    self.options_errorCheck = QCheckBox(self.optionsPage,"options_errorCheck")
    self.options_errorCheck.setChecked(1)
    Layout19.addWidget(self.options_errorCheck)
    Spacer23 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout19.addItem(Spacer23)
    optionsPageLayout.addLayout(Layout19)

    Layout20 = QHBoxLayout(None,0,6,"Layout20")

    self.options_rejectionCheck = QCheckBox(self.optionsPage,"options_rejectionCheck")
    Layout20.addWidget(self.options_rejectionCheck)
    Spacer24 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20.addItem(Spacer24)
    optionsPageLayout.addLayout(Layout20)

    Layout20_2 = QHBoxLayout(None,0,6,"Layout20_2")

    self.options_oobCheck = QCheckBox(self.optionsPage,"options_oobCheck")
    Layout20_2.addWidget(self.options_oobCheck)
    Spacer24_2 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout20_2.addItem(Spacer24_2)
    optionsPageLayout.addLayout(Layout20_2)

    Layout15 = QHBoxLayout(None,0,6,"Layout15")

    self.ButtonGroup1 = QButtonGroup(self.optionsPage,"ButtonGroup1")
    self.ButtonGroup1.setColumnLayout(0,Qt.Vertical)
    self.ButtonGroup1.layout().setSpacing(6)
    self.ButtonGroup1.layout().setMargin(11)
    ButtonGroup1Layout = QVBoxLayout(self.ButtonGroup1.layout())
    ButtonGroup1Layout.setAlignment(Qt.AlignTop)

    self.options_sourceAll = QRadioButton(self.ButtonGroup1,"options_sourceAll")
    self.options_sourceAll.setChecked(1)
    ButtonGroup1Layout.addWidget(self.options_sourceAll)

    self.options_sourceTraining = QRadioButton(self.ButtonGroup1,"options_sourceTraining")
    ButtonGroup1Layout.addWidget(self.options_sourceTraining)

    self.options_sourceHoldout = QRadioButton(self.ButtonGroup1,"options_sourceHoldout")
    ButtonGroup1Layout.addWidget(self.options_sourceHoldout)
    Layout15.addWidget(self.ButtonGroup1)
    Spacer19 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout15.addItem(Spacer19)
    optionsPageLayout.addLayout(Layout15)
    Spacer22 = QSpacerItem(16,16,QSizePolicy.Minimum,QSizePolicy.Expanding)
    optionsPageLayout.addItem(Spacer22)
    self.tabWidget.insertTab(self.optionsPage,QString.fromLatin1(""))

    self.dataPage = QWidget(self.tabWidget,"dataPage")
    self.tabWidget.insertTab(self.dataPage,QString.fromLatin1(""))

    self.filePage = QWidget(self.tabWidget,"filePage")
    self.tabWidget.insertTab(self.filePage,QString.fromLatin1(""))
    ScreenCompositeDialogLayout.addWidget(self.tabWidget)

    Layout46 = QHBoxLayout(None,0,6,"Layout46")

    self.TextLabel3 = QLabel(self,"TextLabel3")
    self.TextLabel3.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Preferred,0,0,self.TextLabel3.sizePolicy().hasHeightForWidth()))
    self.TextLabel3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    Layout46.addWidget(self.TextLabel3)

    self.options_modelName = QLineEdit(self,"options_modelName")
    self.options_modelName.setFrameShadow(QLineEdit.Plain)
    self.options_modelName.setReadOnly(1)
    Layout46.addWidget(self.options_modelName)
    ScreenCompositeDialogLayout.addLayout(Layout46)

    Layout45 = QHBoxLayout(None,0,6,"Layout45")

    self.buttonHelp = QPushButton(self,"buttonHelp")
    self.buttonHelp.setAutoDefault(1)
    Layout45.addWidget(self.buttonHelp)
    Spacer12 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout45.addItem(Spacer12)

    self.buttonRefresh = QPushButton(self,"buttonRefresh")
    self.buttonRefresh.setEnabled(1)
    Layout45.addWidget(self.buttonRefresh)
    Spacer13 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout45.addItem(Spacer13)

    self.buttonScreen = QPushButton(self,"buttonScreen")
    self.buttonScreen.setEnabled(0)
    Layout45.addWidget(self.buttonScreen)
    Spacer13_2 = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout45.addItem(Spacer13_2)

    self.buttonScreenAll = QPushButton(self,"buttonScreenAll")
    self.buttonScreenAll.setEnabled(0)
    Layout45.addWidget(self.buttonScreenAll)
    Spacer15_2 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout45.addItem(Spacer15_2)

    self.buttonClose = QPushButton(self,"buttonClose")
    self.buttonClose.setAutoDefault(1)
    self.buttonClose.setDefault(1)
    Layout45.addWidget(self.buttonClose)
    ScreenCompositeDialogLayout.addLayout(Layout45)

    self.languageChange()

    self.resize(QSize(472,302).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.buttonScreen,SIGNAL("clicked()"),self.screenClick)
    self.connect(self.buttonClose,SIGNAL("clicked()"),self.accept)
    self.connect(self.buttonRefresh,SIGNAL("pressed()"),self.refreshClick)
    self.connect(self.buttonScreenAll,SIGNAL("clicked()"),self.screenAllClick)


  def languageChange(self):
    self.setCaption(self.__tr("ScreenComposite"))
    self.options_thresholdEdit.setText(self.__tr("0.0"))
    QToolTip.add(self.options_thresholdEdit,self.__tr("Screening threshold.  Enables high-confidence predictions.  This should be between 0 and 1."))
    self.TextLabel1.setText(self.__tr("Threshold"))
    self.options_errorCheck.setText(self.__tr("Show Errors"))
    QToolTip.add(self.options_errorCheck,self.__tr("Toggles creation of a worksheet with bad predictions."))
    self.options_rejectionCheck.setText(self.__tr("Show Rejections"))
    self.options_oobCheck.setText(self.__tr("Generalization Error Estimate"))
    QToolTip.add(self.options_oobCheck,self.__tr("Generate the \"Out-of-Bag\" error estimate instead of the normal estimate. This only makes sense when screening the original data set."))
    self.ButtonGroup1.setTitle(self.__tr("Pool"))
    self.options_sourceAll.setText(self.__tr("All Data"))
    QToolTip.add(self.options_sourceAll,self.__tr("Entire data set"))
    self.options_sourceTraining.setText(self.__tr("Training Set"))
    QToolTip.add(self.options_sourceTraining,self.__tr("Only the data used to build the model (only meaningful when screening the same data set used for model construction)."))
    self.options_sourceHoldout.setText(self.__tr("Holdout Set"))
    QToolTip.add(self.options_sourceHoldout,self.__tr("Only the data in the holdout set (only meaningful when screening the same data set used for model construction)."))
    self.tabWidget.changeTab(self.optionsPage,self.__tr("Options"))
    self.tabWidget.changeTab(self.dataPage,self.__tr("Database"))
    self.tabWidget.changeTab(self.filePage,self.__tr("File"))
    self.TextLabel3.setText(self.__tr("Model:"))
    self.buttonHelp.setText(self.__tr("Help"))
    self.buttonHelp.setAccel(self.__tr("F1"))
    self.buttonRefresh.setText(self.__tr("Refresh"))
    QToolTip.add(self.buttonRefresh,self.__tr("Grabs a composite model from our parent."))
    self.buttonScreen.setText(self.__tr("Screen"))
    QToolTip.add(self.buttonScreen,self.__tr("Screen the current composite model."))
    self.buttonScreenAll.setText(self.__tr("Screen All"))
    QToolTip.add(self.buttonScreenAll,self.__tr("Screen all composite models in the CompositeInteract window."))
    self.buttonClose.setText(self.__tr("Close"))
    self.buttonClose.setAccel(QString.null)
    QToolTip.add(self.buttonClose,self.__tr("Close this window."))


  def refreshClick(self):
    print "ScreenCompositeDialog.refreshClick(): Not implemented yet"

  def screenClick(self):
    print "ScreenCompositeDialog.screenClick(): Not implemented yet"

  def screenAllClick(self):
    print "ScreenCompositeDialog.screenAllClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("ScreenCompositeDialog",s,c)
