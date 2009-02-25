#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" provides a window with an interactive python shell

"""

import sys
from qt import *
from pyRDKit.qtGui.forms.pyshell import PyShell as _Form
from code import InteractiveInterpreter

class PyShell(_Form,InteractiveInterpreter):
  """  This class provides a window for interacting with a python intepreter

  **Form:**   _forms.pyshell.PyShell_

  The window is has 2 panes in a splitter:

     1) allows input,

     2) displays output


  The shell maintains a history of commands entered, so doing a log
  would be pretty easy (but that's not implemented yet)
  
  """
  def __init__(self,parent=None,name=None,locals=None):
    """ constructor

    **Arguments**

      - parent: (optional) the parent widget for this one

      - name: (optional) the name of the widget to use

      - locals: (optional) if provided, this should be a dictionary
        which will be passed on to set the namespace for the
        _InteractiveInterpreter_ we launch.

    """
    _Form.__init__(self,parent,name)
    InteractiveInterpreter.__init__(self,locals)
    self.cmds = []
    self._activeCmd = -1
    
  def write(self,txt):
    """ dumps output to the output window

    """
    #self.output.setBold(1)
    self.output.insertParagraph(txt,-1)
    #self.output.setBold(0)

  def doCommand(self,verbose=0):
    """ executes a command

    **Notes**:

       - the command is inserted into the output pane and added to the
         history list.

    """
    txt = str(self.input.text())
    txt = txt.replace('\n','')
    if verbose: print 'cmd:',repr(txt)
    stdout = sys.stdout
    naf = NotAFile(self.output)
    sys.stdout = naf
    if len(txt.replace('\n','')):
      self.cmds.append(txt)
      self._activeCmd = len(self.cmds)
      self.input.clear()
      self.output.setBold(1)
      self.output.insertParagraph(txt+'\n',-1)
      self.output.setBold(0)

    res = self.runsource(txt)
    sys.stdout = stdout
    if verbose: print 'res:',repr(res)
    if not res:
      naf._finish()

  def _updateFromHistory(self):
    """ INTERNAL USE ONLY

     updates the contents of the input pane based upon the current
     history item

    """
    try:
      txt = self.cmds[self._activeCmd]
    except IndexError:
      txt = ''
    self.input.setText(txt)
    self.input.setCursorPosition(0,len(txt))
    
  def ctrlP(self):
    """ callback for Ctrl-P presses: scrolls up in the history

    """
    if self._activeCmd > 0:
      self._activeCmd -= 1
      self._updateFromHistory()

  def ctrlN(self):
    """ callback for Ctrl-N presses: scrolls down in the history

    """
    nCmds = len(self.cmds)
    if self._activeCmd < nCmds:
      self._activeCmd += 1
      self._updateFromHistory()

  def fileExit(self):
    self.hide()
      
class NotAFile:
  """ class used by the _InteractiveInterpreter_ component of the
      _PyShell_ to fake a file

  """
  def __init__(self,where,color=None,font=None):
    self.where = where
    self.color = color
    self.font = font
    self._lines = []
  def write(self,what):
    self._lines.append(what)
  def writelines(self,what):
    self._lines += what
  def flush(self):
    pass
  def _finish(self):
    self.where.insertParagraph(''.join(self._lines),-1)

if __name__ == '__main__':
  import sys
  a = QApplication(sys.argv)
  widg = PyShell()
  widg.show()
  a.setMainWidget(widg)
  a.exec_loop()


