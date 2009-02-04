#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
"""  dialog for establishing db connections

"""    
from pyRDKit import RDConfig
from qt import *
from pyRDKit.qtGui.DbConnWidgetImpl import insertConnWidget
from pyRDKit.qtGui.DbDialog import DbDialog


class DbConnDialog(DbDialog):
  """  dialog for establishing db connections

    This is just a _DbDialog_ containing a _DbConnWidget_
    
  """    
  def __init__(self,parent=None,initDir=''):
    DbDialog.__init__(self,parent,initDir)
    self.setDbWidget(insertConnWidget(self))

