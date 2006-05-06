# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DbQueryWidget.ui'
#
# Created: Mon Jan 23 08:34:47 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class DbQueryWidget(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
    QWidget.__init__(self,parent,name,fl)

    if not name:
      self.setName("DbQueryWidget")


    DbQueryWidgetLayout = QGridLayout(self,1,1,4,2,"DbQueryWidgetLayout")

    self.db_tabelLabel = QLabel(self,"db_tabelLabel")
    self.db_tabelLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_tabelLabel,3,0)

    self.userField = QLineEdit(self,"userField")
    self.userField.setFrameShape(QLineEdit.LineEditPanel)
    self.userField.setFrameShadow(QLineEdit.Sunken)

    DbQueryWidgetLayout.addMultiCellWidget(self.userField,1,1,1,3)

    self.db_sqlLabel_3 = QLabel(self,"db_sqlLabel_3")
    self.db_sqlLabel_3.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sqlLabel_3,1,0)

    self.db_sqlLabel_2_2 = QLabel(self,"db_sqlLabel_2_2")
    self.db_sqlLabel_2_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sqlLabel_2_2,6,0)

    self.db_sqlLabel = QLabel(self,"db_sqlLabel")
    self.db_sqlLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sqlLabel,4,0)

    self.db_sourceLabel = QLabel(self,"db_sourceLabel")
    self.db_sourceLabel.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sourceLabel,0,0)

    self.tableCombo = QComboBox(0,self,"tableCombo")
    self.tableCombo.setEnabled(0)
    self.tableCombo.setAutoCompletion(1)
    self.tableCombo.setDuplicatesEnabled(0)

    DbQueryWidgetLayout.addMultiCellWidget(self.tableCombo,3,3,1,3)

    self.passwordField = QLineEdit(self,"passwordField")
    self.passwordField.setEchoMode(QLineEdit.Password)

    DbQueryWidgetLayout.addMultiCellWidget(self.passwordField,2,2,1,3)

    self.db_sqlLabel_2 = QLabel(self,"db_sqlLabel_2")
    self.db_sqlLabel_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sqlLabel_2,5,0)

    self.db_sqlLabel_3_2 = QLabel(self,"db_sqlLabel_3_2")
    self.db_sqlLabel_3_2.setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)

    DbQueryWidgetLayout.addWidget(self.db_sqlLabel_3_2,2,0)

    self.sqlWhereBox = QLineEdit(self,"sqlWhereBox")
    self.sqlWhereBox.setEnabled(0)

    DbQueryWidgetLayout.addMultiCellWidget(self.sqlWhereBox,5,5,1,3)

    self.sqlJoinBox = QLineEdit(self,"sqlJoinBox")
    self.sqlJoinBox.setEnabled(0)

    DbQueryWidgetLayout.addMultiCellWidget(self.sqlJoinBox,6,6,1,3)

    self.nameBox = QLineEdit(self,"nameBox")

    DbQueryWidgetLayout.addMultiCellWidget(self.nameBox,0,0,1,2)

    self.nameButton = QToolButton(self,"nameButton")

    DbQueryWidgetLayout.addWidget(self.nameButton,0,3)

    self.sqlWhatBox = QLineEdit(self,"sqlWhatBox")
    self.sqlWhatBox.setEnabled(0)

    DbQueryWidgetLayout.addWidget(self.sqlWhatBox,4,1)

    self.sqlChooseButton = QToolButton(self,"sqlChooseButton")
    self.sqlChooseButton.setEnabled(0)

    DbQueryWidgetLayout.addMultiCellWidget(self.sqlChooseButton,4,4,2,3)

    self.languageChange()

    self.resize(QSize(429,169).expandedTo(self.minimumSizeHint()))
    self.clearWState(Qt.WState_Polished)

    self.connect(self.nameButton,SIGNAL("clicked()"),self.dbFileClick)
    self.connect(self.sqlChooseButton,SIGNAL("clicked()"),self.dbChooseClick)


  def languageChange(self):
    self.setCaption(self.__tr("DbQueryWidget"))
    self.db_tabelLabel.setText(self.__tr("Table"))
    self.userField.setText(self.__tr("sysdba"))
    QToolTip.add(self.userField,self.__tr("username for DB access"))
    self.db_sqlLabel_3.setText(self.__tr("User"))
    QToolTip.add(self.db_sqlLabel_3,self.__tr("username for DB access"))
    self.db_sqlLabel_2_2.setText(self.__tr("SQL Join"))
    self.db_sqlLabel.setText(self.__tr("SQL What"))
    self.db_sourceLabel.setText(self.__tr("Database"))
    QToolTip.add(self.tableCombo,self.__tr("name of the database table"))
    self.passwordField.setText(self.__tr("masterkey"))
    QToolTip.add(self.passwordField,self.__tr("password for DB access"))
    self.db_sqlLabel_2.setText(self.__tr("SQL Where"))
    self.db_sqlLabel_3_2.setText(self.__tr("Password"))
    QToolTip.add(self.db_sqlLabel_3_2,self.__tr("password for DB access"))
    QToolTip.add(self.sqlWhereBox,self.__tr("where specification for the SQL query"))
    QToolTip.add(self.sqlJoinBox,self.__tr("where specification for the SQL query"))
    QToolTip.add(self.nameBox,self.__tr("file containing the database"))
    self.nameButton.setText(self.__tr("..."))
    self.sqlWhatBox.setText(self.__tr("*"))
    QToolTip.add(self.sqlWhatBox,self.__tr("which fields are to be taken from the table "))
    self.sqlChooseButton.setText(self.__tr("Choose"))
    QToolTip.add(self.sqlChooseButton,self.__tr("Select the columns from a list."))


  def dbFileClick(self):
    print "DbQueryWidget.dbFileClick(): Not implemented yet"

  def dbChooseClick(self):
    print "DbQueryWidget.dbChooseClick(): Not implemented yet"

  def __tr(self,s,c = None):
    return qApp.translate("DbQueryWidget",s,c)
