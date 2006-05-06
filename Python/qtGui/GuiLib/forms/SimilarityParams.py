# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SimilarityParams.ui'
#
# Created: Mon Jan 23 08:36:05 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class SimilarityParamsWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if not name:
      self.setName("SimilarityParamsWidget")

    self.setEnabled(1)
    self.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.sizePolicy().hasHeightForWidth()))

    SimilarityParamsWidgetLayout = QVBoxLayout(self,4,4,"SimilarityParamsWidgetLayout")

    self.ButtonGroup7 = QButtonGroup(self,"ButtonGroup7")
    self.ButtonGroup7.setEnabled(1)
    self.ButtonGroup7.setColumnLayout(0,Qt.Vertical)
    self.ButtonGroup7.layout().setSpacing(2)
    self.ButtonGroup7.layout().setMargin(4)
    ButtonGroup7Layout = QVBoxLayout(self.ButtonGroup7.layout())
    ButtonGroup7Layout.setAlignment(Qt.AlignTop)

    Layout11 = QHBoxLayout(None,0,6,"Layout11")

    self.fragmentRadio = QRadioButton(self.ButtonGroup7,"fragmentRadio")
    self.fragmentRadio.setChecked(1)
    Layout11.addWidget(self.fragmentRadio)
    Spacer6 = QSpacerItem(80,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout11.addItem(Spacer6)

    self.TextLabel1 = QLabel(self.ButtonGroup7,"TextLabel1")
    self.TextLabel1.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    Layout11.addWidget(self.TextLabel1)

    self.fragmentNumBits = QLineEdit(self.ButtonGroup7,"fragmentNumBits")
    self.fragmentNumBits.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.fragmentNumBits.sizePolicy().hasHeightForWidth()))
    self.fragmentNumBits.setMinimumSize(QSize(60,0))
    self.fragmentNumBits.setMaximumSize(QSize(60,32767))
    Layout11.addWidget(self.fragmentNumBits)

    self.TextLabel1_2_2 = QLabel(self.ButtonGroup7,"TextLabel1_2_2")
    self.TextLabel1_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    Layout11.addWidget(self.TextLabel1_2_2)

    self.fragmentMinSize = QSpinBox(self.ButtonGroup7,"fragmentMinSize")
    self.fragmentMinSize.setMaxValue(10)
    self.fragmentMinSize.setMinValue(1)
    self.fragmentMinSize.setValue(1)
    Layout11.addWidget(self.fragmentMinSize)

    self.TextLabel1_2 = QLabel(self.ButtonGroup7,"TextLabel1_2")
    self.TextLabel1_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
    Layout11.addWidget(self.TextLabel1_2)

    self.fragmentMaxSize = QSpinBox(self.ButtonGroup7,"fragmentMaxSize")
    self.fragmentMaxSize.setMaxValue(10)
    self.fragmentMaxSize.setMinValue(1)
    self.fragmentMaxSize.setValue(7)
    Layout11.addWidget(self.fragmentMaxSize)
    ButtonGroup7Layout.addLayout(Layout11)

    Layout28 = QHBoxLayout(None,0,6,"Layout28")

    self.MACCSRadio = QRadioButton(self.ButtonGroup7,"MACCSRadio")
    self.MACCSRadio.setChecked(0)
    Layout28.addWidget(self.MACCSRadio)
    Spacer8 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout28.addItem(Spacer8)
    ButtonGroup7Layout.addLayout(Layout28)

    Layout25 = QHBoxLayout(None,0,6,"Layout25")

    self.RadioButton19 = QRadioButton(self.ButtonGroup7,"RadioButton19")
    self.RadioButton19.setEnabled(0)
    Layout25.addWidget(self.RadioButton19)
    Spacer7 = QSpacerItem(30,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout25.addItem(Spacer7)
    ButtonGroup7Layout.addLayout(Layout25)

    Layout27 = QGridLayout(None,1,1,0,6,"Layout27")

    self.patternRadio = QRadioButton(self.ButtonGroup7,"patternRadio")
    self.patternRadio.setEnabled(0)

    Layout27.addWidget(self.patternRadio,0,0)

    self.moduleRadio = QRadioButton(self.ButtonGroup7,"moduleRadio")
    self.moduleRadio.setEnabled(0)

    Layout27.addWidget(self.moduleRadio,1,0)

    self.moduleEdit = QLineEdit(self.ButtonGroup7,"moduleEdit")
    self.moduleEdit.setEnabled(0)

    Layout27.addWidget(self.moduleEdit,1,1)

    Layout79_2 = QHBoxLayout(None,0,6,"Layout79_2")

    self.patternEdit = QLineEdit(self.ButtonGroup7,"patternEdit")
    self.patternEdit.setEnabled(0)
    Layout79_2.addWidget(self.patternEdit)

    self.patternFileButton = QToolButton(self.ButtonGroup7,"patternFileButton")
    self.patternFileButton.setEnabled(0)
    Layout79_2.addWidget(self.patternFileButton)

    Layout27.addLayout(Layout79_2,0,1)
    ButtonGroup7Layout.addLayout(Layout27)
    SimilarityParamsWidgetLayout.addWidget(self.ButtonGroup7)

    layout11 = QHBoxLayout(None,0,6,"layout11")

    self.metric = QButtonGroup(self,"metric")
    self.metric.setColumnLayout(0,Qt.Vertical)
    self.metric.layout().setSpacing(2)
    self.metric.layout().setMargin(4)
    metricLayout = QVBoxLayout(self.metric.layout())
    metricLayout.setAlignment(Qt.AlignTop)

    self.euclidianRadio = QRadioButton(self.metric,"euclidianRadio")
    self.euclidianRadio.setEnabled(0)
    self.euclidianRadio.setChecked(0)
    metricLayout.addWidget(self.euclidianRadio)

    self.tanimotoRadio = QRadioButton(self.metric,"tanimotoRadio")
    self.tanimotoRadio.setChecked(1)
    metricLayout.addWidget(self.tanimotoRadio)

    self.diceRadio = QRadioButton(self.metric,"diceRadio")
    metricLayout.addWidget(self.diceRadio)

    self.cosineRadio = QRadioButton(self.metric,"cosineRadio")
    metricLayout.addWidget(self.cosineRadio)

    self.sokalRadio = QRadioButton(self.metric,"sokalRadio")
    metricLayout.addWidget(self.sokalRadio)
    layout11.addWidget(self.metric)

    self.ButtonGroup5 = QButtonGroup(self,"ButtonGroup5")
    self.ButtonGroup5.setColumnLayout(0,Qt.Vertical)
    self.ButtonGroup5.layout().setSpacing(6)
    self.ButtonGroup5.layout().setMargin(11)
    ButtonGroup5Layout = QVBoxLayout(self.ButtonGroup5.layout())
    ButtonGroup5Layout.setAlignment(Qt.AlignTop)

    Layout14 = QGridLayout(None,1,1,0,6,"Layout14")

    self.greaterThanRadio = QRadioButton(self.ButtonGroup5,"greaterThanRadio")

    Layout14.addWidget(self.greaterThanRadio,1,0)

    self.topNRadio = QRadioButton(self.ButtonGroup5,"topNRadio")
    self.topNRadio.setChecked(1)

    Layout14.addWidget(self.topNRadio,0,0)

    self.greaterThan = QLineEdit(self.ButtonGroup5,"greaterThan")
    self.greaterThan.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.greaterThan.sizePolicy().hasHeightForWidth()))
    self.greaterThan.setMinimumSize(QSize(60,0))
    self.greaterThan.setMaximumSize(QSize(60,32767))

    Layout14.addWidget(self.greaterThan,1,1)

    self.topN = QLineEdit(self.ButtonGroup5,"topN")
    self.topN.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.topN.sizePolicy().hasHeightForWidth()))
    self.topN.setMinimumSize(QSize(60,0))
    self.topN.setMaximumSize(QSize(60,32767))

    Layout14.addWidget(self.topN,0,1)
    ButtonGroup5Layout.addLayout(Layout14)
    layout11.addWidget(self.ButtonGroup5)
    Spacer12 = QSpacerItem(104,16,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout11.addItem(Spacer12)
    SimilarityParamsWidgetLayout.addLayout(layout11)

    self.languageChange()

    self.resize(QSize(445,289).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.patternFileButton,SIGNAL("pressed()"),self.patternFileButtonPressed)


  def languageChange(self):
    self.setCaption(self.__tr("Similarity"))
    self.ButtonGroup7.setTitle(self.__tr("Fingerprint"))
    self.fragmentRadio.setText(self.__tr("Fragments"))
    QToolTip.add(self.fragmentRadio,self.__tr("A binary, fragment-based fingerprinting scheme (similar to Daylight fingerprints)."))
    self.TextLabel1.setText(self.__tr("Num Bits"))
    self.fragmentNumBits.setText(self.__tr("2048"))
    QToolTip.add(self.fragmentNumBits,self.__tr("Number of bits to include in the fingerprint"))
    self.TextLabel1_2_2.setText(self.__tr("Min Path"))
    QToolTip.add(self.fragmentMinSize,self.__tr("Minimum fragment length to be included in the fingerprint"))
    self.TextLabel1_2.setText(self.__tr("Max Path"))
    QToolTip.add(self.fragmentMaxSize,self.__tr("Maximum fragment length to be included in the fingerprint"))
    self.MACCSRadio.setText(self.__tr("MACCS Keys"))
    QToolTip.add(self.MACCSRadio,self.__tr("A binary fingerprint based on the published MACCS keys"))
    self.RadioButton19.setText(self.__tr("EState"))
    QToolTip.add(self.RadioButton19,self.__tr("Use Hall-Kier electrotopological state fingerprints."))
    self.patternRadio.setText(self.__tr("Pattern"))
    QToolTip.add(self.patternRadio,self.__tr("Define the fingerprint using a set of SMARTS patterns in a file."))
    self.moduleRadio.setText(self.__tr("Module"))
    QToolTip.add(self.moduleRadio,self.__tr("Use a fingerprint calculation function from a module, specified here as: Module.Function"))
    QToolTip.add(self.moduleEdit,self.__tr("Use a fingerprint calculation function from a module, specified here as: Module.Function"))
    QToolTip.add(self.patternEdit,self.__tr("file containing the SMARTS patterns"))
    self.patternFileButton.setText(self.__tr("..."))
    self.metric.setTitle(self.__tr("Similarity Metric"))
    self.euclidianRadio.setText(self.__tr("Euclidian"))
    QToolTip.add(self.euclidianRadio,self.__tr("Use the Euclidian distance between vectors."))
    self.tanimotoRadio.setText(self.__tr("Tanimoto"))
    QToolTip.add(self.tanimotoRadio,self.__tr("Use the Tanimoto distance between vectors (only makes sense for bit vectors)."))
    self.diceRadio.setText(self.__tr("Dice"))
    QToolTip.add(self.diceRadio,self.__tr("Use the DICE distance between vectors (only makes sense for bit vectors)."))
    self.cosineRadio.setText(self.__tr("Cosine"))
    QToolTip.add(self.cosineRadio,self.__tr("Use the cosine distance between vectors (only makes sense for bit vectors)."))
    self.sokalRadio.setText(self.__tr("Sokal"))
    QToolTip.add(self.sokalRadio,self.__tr("Use the Sokal distance between vectors (only makes sense for bit vectors)."))
    self.ButtonGroup5.setTitle(self.__tr("Results"))
    self.greaterThanRadio.setText(self.__tr("Greater Than"))
    QToolTip.add(self.greaterThanRadio,self.__tr("All compounds with a similarity above the threshold will be retrieved."))
    self.topNRadio.setText(self.__tr("Top N"))
    QToolTip.add(self.topNRadio,self.__tr("Get only the top N most similar hits.  (N=-1 gives a ranked list of all hits)"))
    self.greaterThan.setText(self.__tr("0.75"))
    QToolTip.add(self.greaterThan,self.__tr("All compounds with a similarity above the threshold will be retrieved."))
    self.topN.setText(self.__tr("10"))
    QToolTip.add(self.topN,self.__tr("Get only the top N most similar hits.  (N=-1 gives a ranked list of all hits)"))


  def patternFileButtonPressed(self):
    print "SimilarityParamsWidget.patternFileButtonPressed(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("SimilarityParamsWidget",s,c)
