# $Id$
#----------------------------------------------------------------------
# Name:        Logger.py
# Purpose:     Provides a Logger class which can be wrapped around another
#              python class to log method calls and attribute changes.
# Requires:    Python 2.0 (or higher?)
#
# Author:      greg Landrum (Landrum@RationalDiscovery.com)
# License:
#
#       Copyright (c) 2001-2006 Greg Landrum and Rational Discovery LLC,
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
#----------------------------------------------------------------------
""" Provides a Logger class which can be wrapped around another
             python class to log method calls and attribute changes.
"""

import types
import re

reType = type(re.compile('f*'))
stringTypes = [types.StringType,types.UnicodeType,reType]
               
def isPresent(what,checkList):
  """ checks to see if any of the regular expressions in a list match a string

   **Arguments**

     - what: the thing to match against

     - checkList: the list of regexps to try the match with

   **Returns**

     1 or 0 depending upon whether or not _what_ was found in checkList

   **Notes**

     - the search is done using _match_, so we match *any* part of _what_

  """
  for entry in checkList:
    if type(entry) == reType:
      if entry.match(what) is not None:
        return 1
    else:
      if what == entry:
        return 1
  return 0

class Callable(object):
  """ this is just a simple class with a __call__ method we'll use to pass
      any method invocations back to the caller.

  """
  def __init__(self,log,obj,method):
    """ Constructor
    
      **Arguments:**

        - log: a python list (or anything else supporting the append method)
             which is used to log the actual invocation of a method

        - obj: the object which is calling the method

        - method: the *name* of the method to be invoked

    """
    self._obj = obj
    self._log = log
    self._method = method
  def __call__(self,*args,**kwargs):
    """ logs the method name and arguments and makes the call

    """
    self._log.append((self._method,args,kwargs))
    return getattr(self._obj,self._method)(*args,**kwargs)

class Logger(object):
  """ This is the actual wrapper class.

    The wrapper is fairly thin; it only has one methods of its own:

      - _LoggerGetLog()

    and then several instance variables:

      - _loggerFlushCommand

      - _loggerClass

      - _loggerObj

      - _loggerCommandLog

      - _loggerIgnore

    These names were chosen to minimize the likelihood of a collision
     with the attributes of a wrapped class.  Obviously... ;-)

    The general idea of using this is that you wrap a class in the logger,
      and then use the class as you normally would.  Whenever you want to 
      get the contents of the log (for example after running your program for
      a while), you can call _loggerCommandLog.  The resulting list can be
      played back in another (or the same) object using the replay() function
      defined below.

    The logger can, optionally, be set to flush its log whenever a method with
      a particular name is invoked.  For example, you may want to be wrapping
      some kind of drawing canvas and want to reset the log whenever the canvas
      is cleared because there's no point in storing commands which will have
      no effect on the final drawing.

    **Note**

      because of the way I've worked this, the log will actually be flushed
      whenever the client program accesses the flush method, it need not be invoked.
      i.e. if the loggerFlushCommand is 'foo', then doing either wrappedObj.foo() or
      wrappedObj.foo will reset the log.  This is undesirable and will most likely
      be fixed in a future version

  """
  def __init__(self,klass,*args,**kwargs):
    """  Constructor

      **Arguments**
      
        The one required argument here is _klass_, which is the class
        to be wrapped.

      **Optional Keyword Arguments**

        - loggerFlushCommand: the name of the attribute which will flush the log

        - loggerIgnore: a list of regexps defining methods which should not be
          logged

      **All other arguments are passed to the constructor for _klass_ **
      
    """
    if kwargs.has_key('loggerFlushCommand'):
      self.__dict__['_loggerFlushCommand'] = kwargs['loggerFlushCommand']
      del kwargs['loggerFlushCommand']
    else:
      self.__dict__['_loggerFlushCommand'] = None
    if kwargs.has_key('loggerIgnore'):
      tempL = kwargs['loggerIgnore']
      for entry in tempL:
        if type(entry) not in stringTypes:
          raise ValueError('All entries in loggerIgnore must be either strings or regexps')
      self.__dict__['_loggerIgnore'] = tempL
      del kwargs['loggerIgnore']
    else:
      self.__dict__['_loggerIgnore'] = []
    self.__dict__['_loggerClass'] = klass
    self.__dict__['_loggerObj'] = klass(*args,**kwargs)
    self.__dict__['_loggerCommandLog'] = []

  def _LoggerGetLog(self):
    """ Returns the contents of the command log as a python list

    """
    return self._loggerCommandLog

  def __getattr__(self,which):
    """ here's where the logging of method invocations takes place

    """
    if hasattr(self._loggerObj,which):
      tmpAttr = getattr(self._loggerObj,which)
      if type(tmpAttr) == types.MethodType:
        loggerFlushCommand = self._loggerFlushCommand
        if which == loggerFlushCommand:
          self._loggerCommandLog = []
          return Callable([],self._loggerObj,which)
        elif self._loggerIgnore != [] and isPresent(which,self._loggerIgnore):
          return Callable([],self._loggerObj,which)
        else:
          return Callable(self._loggerCommandLog,self._loggerObj,which)
      else:
        return tmpAttr
    else:
      raise AttributeError('%s instance has no attribute %s'%(repr(self._loggerClass.__name__),repr(which)))

  def __setattr__(self,which,val):
    """ setattr calls (i.e. wrappedObject.foo = 1) are also logged

    """
    d = self.__dict__
    if d.has_key(which):
      d[which] = val
    else:
      self._loggerCommandLog.append((setattr,which,val))
      setattr(self._loggerObj,which,val)
      return val
    
def replay(logItems,obj):
  """ loops through the items in a Logger log list and invokes
     them in obj

    **Arguments**

      - logItems: a list of 3 tuples containing:

         1) method name

         2) tuple of positional arguments

         3) dictionary of keyword arguments

      - obj: the object in which the log should be replayed  

    **Returns**

      a list with the the return values of all the method
      invocations/attribute assignments

  """
  if isinstance(logItems,Logger):
    logItems = logItems._LoggerGetLog()
  resList = []
  for method,a1,a2 in logItems:
    if callable(method):
      if method == setattr:
        method(obj,a1,a2)
        resList.append(a2)
      else:
        resList.append(method(obj,a1,a2))        
    else:
      a = getattr(obj,method)
      resList.append(apply(a,a1,a2))
  return resList 

