#
#  Copyright (C) 2000  greg Landrum
#
""" Everything required to load and initialize mixins

  **What's in a mixin?**
  
    1) REQUIRED_MIXINS: a list of the names (strings) of other mixins this
                     mixin requires.

    2) MODULES_ALTERED: a list of the names (strings) of the modules this
                     mixin needs so that their contents can be modified.

    3) METHODS_DEFINED: a dictionary consisting of 'function name':'final name'
                     pairs where 'function name' is the name of the function
                     in this module and 'final name' is the fully qualified
                     (module and class) name of where it ends up.

    4) VARS_TO_SAVE: a list of variable names which should be saved when the
                     _StateSaver_ is invoked  

    5) LocalInit:       function containing initialization code which is called
                     when after all mixins have been loaded and the final class
                     has been instantiated.
"""
import sys

def _OrderMixIns(mixInDict,loadOrder):
  """ determines the overall loadOrder of the mixins

    **Arguments**

      - mixInDict: a dictionary of mixins to be loaded as keys and their
         requirements as values

      - loadOrder: a list containing the order in which the
        mixins should be loaded

    **Notes**

      - _loadOrder_ is modified during this process
      
  """    
  localDict = mixInDict.copy()
  tempDict = localDict.copy()
  while len(localDict):
    keys = localDict.keys()
    for key in keys:
      if len(localDict[key]) == 0:
        loadOrder.append(key)
        # remove this key from everyone else's dependency lists
        for key2 in keys:
          try:
            tempDict[key2].remove(key)
          except (KeyError, ValueError):
            pass
          
        # remove the zero-length node
        del tempDict[key]
    localDict = tempDict
            
        
def _FindRequiredMixIns(mixInList,mixInDict):
  """Determines the complete list of mixins which must be loaded

    **Arguments**

      - mixInList: a list of mixins to be loaded

      - mixInDict: a dictionary of mixins to be loaded as keys and their
         requirements as values

    **Notes**

      - _mixInDict_ is modified during this process
      
  """
  for mixIn in mixInList:
    try:
      mod = __import__(mixIn)
    except ImportError:
      import traceback
      traceback.print_exc()
      sys.stderr.write('\n\nERROR: Mixin %s not found, exiting\n'%(repr(mixIn)))
      sys.exit(1)

    if hasattr(mod,'REQUIRED_MIXINS'):
      reqs = mod.REQUIRED_MIXINS
      for req in reqs:
        if not mixInDict.has_key(req):
          _FindRequiredMixIns([req], mixInDict)
      mixInDict[mixIn] = reqs
    else:
      print 'mixIn %s has no REQUIRED_MIXINS'%mixIn

def _FindLoadOrder(mixInList,loadOrder=[],mixInDict={}):
  """Determines the complete list of mixins which must be loaded
     and in what order that should happen.

    **Arguments**

      - mixInList: the initial list of mixins to be loaded

      - loadOrder: a list which will be used to return the load order

      - mixInDict: a dictionary of mixins to be loaded as keys and their
         requirements as values

    **Notes**

      - _loadOrder_ is modified during this process

      -we try to be tolerant of missing mixin requirements: error
          messages are printed. This tolerance may be a bad idea.

  """

  # first get the list of all requirements
  mixInDict = {}
  _FindRequiredMixIns(mixInList,mixInDict)

  # now order them
  _OrderMixIns(mixInDict,loadOrder)



def LoadMixIns(mixInList,loadOrder=[],searchDirs=[]):
  """load a list of mixins and everything they require

    **Arguments**

      - mixInList: the initial list of mixins to be loaded

      - loadOrder is modified to contain the
         names of all the mixins loaded and the order in which
         the loading has been done

      - searchDirs contains a list of directory
        names which should be loaded to sys.path to help find the mixins

    **Returns**

      a list containing the modules (returned by __import__)
      is returned.

    **Notes**

          -sys.path IS ALTERED HERE.  _InitMixIns_ sets it back to normal.

          -we try to be tolerant of missing mixin requirements: error
            messages are printed. This tolerance may be a bad idea.

  """
  sys.origPath = sys.path[:]
  sys.path = searchDirs + sys.path

  mD = {}
  _FindLoadOrder(mixInList,loadOrder,mD)
  mixInsLoaded = []
  for mixIn in loadOrder:
    mixInModule = __import__(mixIn)
    mixInsLoaded.append(mixInModule)
    if not hasattr(mixInModule,'MODULES_ALTERED'):
      print 'mixIn %s alters no modules'%mixIn
    else:
      mods = mixInModule.MODULES_ALTERED
      for mod in mods:
        exec('import %s'%mod)
    if not hasattr(mixInModule,'METHODS_DEFINED'):
      print 'mixIn %s has no methods'%mixIn      
    else:
      methods = mixInModule.METHODS_DEFINED
      for method in methods.keys():
        val = methods[method]
        exec('%s = mixInModule.%s'%(val,method))

  return mixInsLoaded

def InitMixIns(obj,mixInModuleList):
  """go through the modules (from __import__) in mixInModuleList
     in order and call their _LocalInit_ functions.

     **Arguments**

       - obj: an instantiated object which will be modified by the mixins

       - mixInModuleList: the list of mixins (in order) which should be loaded

     **Notes**  

        this also undoes the modifications to sys.path which
           were made in _LoadMixIns_

  """
  for mixInModule in mixInModuleList:
    if not hasattr(mixInModule,'LocalInit'):
      print 'mixIn %s has no LocalInit()'%mixInModule.__name__
    else:
      init = mixInModule.LocalInit
      init(obj)
  obj._mixInList = mixInModuleList
  sys.path= sys.origPath[:]


if __name__ == '__main__':
  loadOrder = []
  _FindLoadOrder(['teste'],loadOrder)
  print 'loadOrder: ', loadOrder
  
