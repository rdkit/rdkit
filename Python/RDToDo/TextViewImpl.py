# $Id: TextViewImpl.py 4899 2006-01-12 16:20:06Z glandrum $
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time
from forms.textview import TextView as _Form

class TextView(_Form):
  def __init__(self,src,
               parent=None,fl=0):
    _Form.__init__(self,parent=parent,fl=fl)
    self.textBrowser.setText(src)
    self.textBrowser.setReadOnly(1)

