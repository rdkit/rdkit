#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for DbConnWidgets

"""    
from rdkit import RDConfig
from qt import *
from rdkit.qtGui.forms.DbConnWidget import DbConnWidget as _Form
from rdkit.qtGui import DbWidget
import os,os.path
from rdkit.Dbase import DbConnection


def insertConnWidget(parent,*args,**kwargs):
  """ constructs a _DbConnWidget_ and inserts it into a parent

   This uses _DbWidget.insertDbWidget_ to do its work
   
  """
  return DbWidget.insertDbWidget(parent,DbConnWidget,*args,**kwargs)

class DbConnWidget(DbWidget.DbWidgetMixin,_Form):
  """  A simple widget for specifying parameters to set up a
   connection to a database

   The widget includes fields for:

     - the database

     - username

     - password

     - name of the table

   **Form:**  _DbConnWidget_

   **Notes**

     - most of the real functionality here comes from the
       _DbWidget.DbWidgetMixin_ mixin class

  """
  def __init__(self,parent=None,initDir='',clickCallback=None):
    _Form.__init__(self,parent)
    DbWidget.DbWidgetMixin.__init__(self,parent,initDir,clickCallback)

