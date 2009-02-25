#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" contains simple dialog box for holding db widgets

"""    
from pyRDKit import RDConfig
from qt import *

class DbDialog(QDialog):
  """ your basic 3 button dialog box with space for a widget

  """
  def __init__(self,parent=None,initDir=''):
    QDialog.__init__(self,parent)
    self._widget = None
    layout = QVBoxLayout(self)
    
    layoutWidget = QWidget(self,"layout1")
    Layout1 = QHBoxLayout(layoutWidget,0,6,"Layout1")
    self.buttonHelp = QPushButton(layoutWidget,"buttonHelp")
    self.buttonHelp.setText(self.trUtf8("Help"))
    self.buttonHelp.setAccel(4144)
    self.buttonHelp.setAutoDefault(1)
    Layout1.addWidget(self.buttonHelp)
    spacer = QSpacerItem(20,0,QSizePolicy.Expanding,QSizePolicy.Minimum)
    Layout1.addItem(spacer)

    self.buttonOk = QPushButton(layoutWidget,"buttonOk")
    self.buttonOk.setText(self.trUtf8("OK"))
    self.buttonOk.setAccel(0)
    self.buttonOk.setAutoDefault(1)
    self.buttonOk.setDefault(1)
    Layout1.addWidget(self.buttonOk)

    self.buttonCancel = QPushButton(layoutWidget,"buttonCancel")
    self.buttonCancel.setText(self.trUtf8("Cancel"))
    self.buttonCancel.setAccel(0)
    self.buttonCancel.setAutoDefault(1)
    Layout1.addWidget(self.buttonCancel)

    layout.addWidget(layoutWidget)

    self.connect(self.buttonOk,SIGNAL("clicked()"),self,SLOT("accept()"))
    self.connect(self.buttonCancel,SIGNAL("clicked()"),self,SLOT("reject()"))
    self.connect(self.buttonHelp,SIGNAL("clicked()"),self.help)


  def setDbWidget(self,widget):
    """ setter for dbWidget attribute """
    self._widget = widget
    
  def dbWidget(self):
    """ returns our widget

    """
    return self._widget
  
  def help(self):
    print 'not yet implemented'
    

