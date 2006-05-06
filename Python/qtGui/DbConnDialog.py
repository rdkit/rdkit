#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
"""  dialog for establishing db connections

"""    
import RDConfig
from qt import *
from qtGui.DbConnWidgetImpl import insertConnWidget
from qtGui.DbDialog import DbDialog


class DbConnDialog(DbDialog):
  """  dialog for establishing db connections

    This is just a _DbDialog_ containing a _DbConnWidget_
    
  """    
  def __init__(self,parent=None,initDir=''):
    DbDialog.__init__(self,parent,initDir)
    self.setDbWidget(insertConnWidget(self))

