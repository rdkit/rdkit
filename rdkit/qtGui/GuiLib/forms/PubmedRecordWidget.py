# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'PubmedRecordWidget.ui'
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

class PubmedRecord(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    self.image0 = QPixmap(image0_data)

    if not name:
      self.setName("PubmedRecord")

    self.setIcon(self.image0)

    PubmedRecordLayout = QVBoxLayout(self,2,6,"PubmedRecordLayout")

    layout3 = QGridLayout(None,1,1,0,6,"layout3")

    self.textLabel4_2_4 = QLabel(self,"textLabel4_2_4")
    self.textLabel4_2_4.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout3.addWidget(self.textLabel4_2_4,2,0)

    self.sourceBox = QLineEdit(self,"sourceBox")
    self.sourceBox.setReadOnly(1)

    layout3.addWidget(self.sourceBox,2,1)

    self.authorBox = QLineEdit(self,"authorBox")
    self.authorBox.setReadOnly(1)

    layout3.addWidget(self.authorBox,1,1)

    self.textLabel4_2_2 = QLabel(self,"textLabel4_2_2")
    self.textLabel4_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout3.addWidget(self.textLabel4_2_2,0,0)

    self.titleBox = QLineEdit(self,"titleBox")
    self.titleBox.setReadOnly(1)

    layout3.addWidget(self.titleBox,0,1)

    self.textLabel4_2 = QLabel(self,"textLabel4_2")
    self.textLabel4_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout3.addWidget(self.textLabel4_2,1,0)
    PubmedRecordLayout.addLayout(layout3)

    layout9 = QGridLayout(None,1,1,0,6,"layout9")

    self.textLabel4_2_3_2 = QLabel(self,"textLabel4_2_3_2")
    self.textLabel4_2_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout9.addWidget(self.textLabel4_2_3_2,0,2)

    self.textLabel4_2_3_3 = QLabel(self,"textLabel4_2_3_3")
    self.textLabel4_2_3_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout9.addWidget(self.textLabel4_2_3_3,1,0)

    self.pageBox = QLineEdit(self,"pageBox")
    self.pageBox.setReadOnly(1)

    layout9.addWidget(self.pageBox,0,3)

    self.textLabel4_2_3 = QLabel(self,"textLabel4_2_3")
    self.textLabel4_2_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout9.addWidget(self.textLabel4_2_3,0,0)

    self.idBox = QLineEdit(self,"idBox")
    self.idBox.setReadOnly(1)

    layout9.addWidget(self.idBox,1,3)

    self.textLabel4 = QLabel(self,"textLabel4")
    self.textLabel4.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    layout9.addWidget(self.textLabel4,1,2)

    self.volumeBox = QLineEdit(self,"volumeBox")
    self.volumeBox.setReadOnly(1)

    layout9.addWidget(self.volumeBox,0,1)

    self.yearBox = QLineEdit(self,"yearBox")
    self.yearBox.setReadOnly(1)

    layout9.addWidget(self.yearBox,1,1)
    PubmedRecordLayout.addLayout(layout9)

    self.abstractBox = QGroupBox(self,"abstractBox")
    self.abstractBox.setColumnLayout(0,Qt.Vertical)
    self.abstractBox.layout().setSpacing(2)
    self.abstractBox.layout().setMargin(2)
    abstractBoxLayout = QVBoxLayout(self.abstractBox.layout())
    abstractBoxLayout.setAlignment(Qt.AlignTop)

    self.abstractViewer = QTextBrowser(self.abstractBox,"abstractViewer")
    self.abstractViewer.setTextFormat(QTextBrowser.AutoText)
    abstractBoxLayout.addWidget(self.abstractViewer)
    PubmedRecordLayout.addWidget(self.abstractBox)

    layout7 = QHBoxLayout(None,0,6,"layout7")

    self.keywordBox = QGroupBox(self,"keywordBox")
    self.keywordBox.setColumnLayout(0,Qt.Vertical)
    self.keywordBox.layout().setSpacing(2)
    self.keywordBox.layout().setMargin(2)
    keywordBoxLayout = QVBoxLayout(self.keywordBox.layout())
    keywordBoxLayout.setAlignment(Qt.AlignTop)

    self.keywordViewer = QTextBrowser(self.keywordBox,"keywordViewer")
    keywordBoxLayout.addWidget(self.keywordViewer)
    layout7.addWidget(self.keywordBox)

    self.chemicalBox = QGroupBox(self,"chemicalBox")
    self.chemicalBox.setLineWidth(1)
    self.chemicalBox.setColumnLayout(0,Qt.Vertical)
    self.chemicalBox.layout().setSpacing(2)
    self.chemicalBox.layout().setMargin(2)
    chemicalBoxLayout = QVBoxLayout(self.chemicalBox.layout())
    chemicalBoxLayout.setAlignment(Qt.AlignTop)

    self.chemicalViewer = QTextBrowser(self.chemicalBox,"chemicalViewer")
    chemicalBoxLayout.addWidget(self.chemicalViewer)
    layout7.addWidget(self.chemicalBox)
    PubmedRecordLayout.addLayout(layout7)

    layout8 = QHBoxLayout(None,0,6,"layout8")

    self.saveButton = QPushButton(self,"saveButton")
    layout8.addWidget(self.saveButton)
    spacer4 = QSpacerItem(40,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout8.addItem(spacer4)
    PubmedRecordLayout.addLayout(layout8)

    self.languageChange()

    self.resize(QSize(544,521).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.saveButton,SIGNAL("clicked()"),self.saveButtonClicked)


  def languageChange(self):
    self.setCaption(self.__tr("Pubmed Record"))
    self.textLabel4_2_4.setText(self.__tr("Source"))
    self.textLabel4_2_2.setText(self.__tr("Title"))
    self.textLabel4_2.setText(self.__tr("Authors"))
    self.textLabel4_2_3_2.setText(self.__tr("Pages"))
    self.textLabel4_2_3_3.setText(self.__tr("Year"))
    self.textLabel4_2_3.setText(self.__tr("Volume"))
    self.textLabel4.setText(self.__tr("PubMedId"))
    self.abstractBox.setTitle(self.__tr("Abstract"))
    self.keywordBox.setTitle(self.__tr("MeSH Keywords"))
    self.chemicalBox.setTitle(self.__tr("MeSH Chemicals"))
    self.saveButton.setText(self.__tr("Save"))
    QToolTip.add(self.saveButton,self.__tr("Write this record to a text (XML) file."))


  def saveButtonClicked(self):
    print "PubmedRecord.saveButtonClicked(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("PubmedRecord",s,c)
