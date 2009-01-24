#
#  Copyright (C) 2000  greg Landrum
#
""" Allows mixins to save their state

  This is done by pickling the variables named in each mixin's *VARS_TO_SAVE*
  variable and writing them out into a file with stub logic for reloading
  the mixins and their variables

"""
import cPickle
import string

_STATE_SAVER_HEADER="""
#
#  Copyright (C) 2001  greg Landrum
#
"""

_STATE_SAVER_STUB="""
def RestoreState(self,varDictString):
  import cPickle
  varDict = cPickle.loads(varDictString)
  for entry in varDict.keys():
    comm = '%s = varDict["%s"]'%(entry,entry)
    exec(comm)
    

if __name__ == '__main__':
  from pyRDKit import RDConfig
  from pyRDKit.GuiFramework import GuiBase

  from pyRDKit.mixins import Loader
  loadOrder = []
  loadedMixIns=Loader.LoadMixIns(REQUIRED_MIXINS,loadOrder,
                                 [RDConfig.RDCodeDir+'/GuiFramework/GuiLib'])
  app = GuiBase.BaseGuiApp(0)
  Loader.InitMixIns(app.frame,loadedMixIns)
  app.frame.FinishMenuBar()
  RestoreState(app.frame,VAR_DICT)
  app.MainLoop()
"""

def _CollectVars(mixInsList,self):
  """ Collects the overall list of *VARS_TO_SAVE*

    The variables are collected into a dictionary which is pickled into a
    string.  The resulting string is returned.

    **Arguments**

      - mixInsList: a list of loaded mixins

      - self: an instantiate object which provides the context for
        'self' references in *VARS_TO_SAVE*

    **Returns**

      a string containing the pickled variable dictionary

    **Notes**

      - if any of the *VARS_TO_SAVE* cannot be pickled, this will generate an error.
      
  """
  stateDict={}
  for module in mixInsList:
    try:
      varsToSave = module.VARS_TO_SAVE
    except:
      pass
    else:
      for var in varsToSave:
        v = eval(var)
        stateDict[var] = v
  res = cPickle.dumps(stateDict)
  return res

def SaveState(mixInsList,self):
  """ saves the state of an object and returns the resulting save file

    **Arguments**

      - mixInsList: a list of loaded mixins

      - self: an instantiate object which provides the context for
        'self' references in *VARS_TO_SAVE*

    **Returns**

      a string containing the contents of the save file

  """
  lines = [_STATE_SAVER_HEADER]
  varDict = _CollectVars(mixInsList,self)
  lines.append('VAR_DICT=%s'%(repr(varDict)))
  mixInNames = map(lambda x: x.__name__, self._mixInList)
  lines.append('REQUIRED_MIXINS=%s'%(repr(mixInNames)))
  lines.append(_STATE_SAVER_STUB)
  return string.join(lines,'\n')
