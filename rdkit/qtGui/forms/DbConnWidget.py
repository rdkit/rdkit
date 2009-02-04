# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DbConnWidget.ui'
#
# Created: Mon Jan 23 08:34:47 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class DbConnWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if not name:
      self.setName("DbConnWidget")


    DbConnWidgetLayout = QGridLayout(self,1,1,2,4,"DbConnWidgetLayout")

    self.db_tabelLabel = QLabel(self,"db_tabelLabel")
    self.db_tabelLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbConnWidgetLayout.addWidget(self.db_tabelLabel,3,0)

    self.userField = QLineEdit(self,"userField")
    self.userField.setFrameShape(QLineEdit.LineEditPanel)
    self.userField.setFrameShadow(QLineEdit.Sunken)

    DbConnWidgetLayout.addMultiCellWidget(self.userField,1,1,1,2)

    self.db_sqlLabel_3 = QLabel(self,"db_sqlLabel_3")
    self.db_sqlLabel_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbConnWidgetLayout.addWidget(self.db_sqlLabel_3,1,0)

    self.db_tabelLabel_2 = QLabel(self,"db_tabelLabel_2")
    self.db_tabelLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbConnWidgetLayout.addWidget(self.db_tabelLabel_2,0,0)

    self.tableCombo = QComboBox(0,self,"tableCombo")
    self.tableCombo.setEnabled(0)
    self.tableCombo.setAutoCompletion(1)
    self.tableCombo.setDuplicatesEnabled(0)

    DbConnWidgetLayout.addMultiCellWidget(self.tableCombo,3,3,1,2)

    self.passwordField = QLineEdit(self,"passwordField")
    self.passwordField.setEchoMode(QLineEdit.Password)

    DbConnWidgetLayout.addMultiCellWidget(self.passwordField,2,2,1,2)

    self.db_sqlLabel_3_2 = QLabel(self,"db_sqlLabel_3_2")
    self.db_sqlLabel_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbConnWidgetLayout.addWidget(self.db_sqlLabel_3_2,2,0)

    self.nameBox = QLineEdit(self,"nameBox")

    DbConnWidgetLayout.addWidget(self.nameBox,0,1)

    self.nameButton = QToolButton(self,"nameButton")

    DbConnWidgetLayout.addWidget(self.nameButton,0,2)

    self.languageChange()

    self.resize(QSize(261,105).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.nameButton,SIGNAL("clicked()"),self.dbFileClick)


  def languageChange(self):
    self.setCaption(self.__tr("DbConnWidget"))
    self.db_tabelLabel.setText(self.__tr("Table"))
    self.userField.setText(self.__tr("sysdba"))
    QToolTip.add(self.userField,self.__tr("username for DB access"))
    self.db_sqlLabel_3.setText(self.__tr("User"))
    QToolTip.add(self.db_sqlLabel_3,self.__tr("username for DB access"))
    self.db_tabelLabel_2.setText(self.__tr("Database"))
    QToolTip.add(self.tableCombo,self.__tr("name of the database table"))
    self.passwordField.setText(self.__tr("masterkey"))
    QToolTip.add(self.passwordField,self.__tr("password for DB access"))
    self.db_sqlLabel_3_2.setText(self.__tr("Password"))
    QToolTip.add(self.db_sqlLabel_3_2,self.__tr("password for DB access"))
    QToolTip.add(self.nameBox,self.__tr("file containing the database"))
    self.nameButton.setText(self.__tr("..."))


  def dbFileClick(self):
    print "DbConnWidget.dbFileClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("DbConnWidget",s,c)
