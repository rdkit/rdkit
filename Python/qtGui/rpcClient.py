#
#  Copyright (C) 2002-2004  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" components required to add an xml-rpc server to a qt application

"""
import sys,exceptions
from qt import *
_rpcEventID=QEvent.User + 2342
_rpcPENDING = 0
_rpcDONE = 1
_rpcEXCEPT = 2

class NoServer(exceptions.RuntimeError):
  pass
class RPCError(exceptions.RuntimeError):
  pass


def HandleRPCEvent(obj,evt,verbose=1):
  """ dispatch function
  
  **Arguments**

    - obj: the object to which the rpc call should be applied

    - evt: the event to use. If this is not an _rpcEvent_ we return.

    - verbose: (optional) selects the verbosity level of the processing


  **Notes**

    - we resolve method names using successive calls to _getattr()_
    
    - there's not much security here
    
    - if no relevant method can be found, we throw an exception

    - information about exceptions thrown by the called method is
      returned to the client to the extent that is possible.

  """
  if evt.type() != _rpcEventID:
    return
  
  if verbose:
    print 'HandleRPCEvent:'
    print '\tmethod:',evt.method()
    print '\tparams:',evt.params()


  # grab the status lock
  evt.rpcStatusLock.lock()
  method = evt.method()
  params = evt.params()
  try:
    # find the actual method by repeatedly applying getattr
    methsplit = method.split('.')
    meth = obj
    for piece in methsplit:
      meth = getattr(meth,piece)
  except AttributeError,msg:
    evt.result = 'No Such Method',msg
    evt.rpcStatus = _rpcEXCEPT
  else:
    try:
      if verbose: print 'apply'
      res = apply(meth,params)
    except:
      import sys,traceback
      traceback.print_exc()
      evt.result = sys.exc_info()[:2]
      evt.rpcStatus = _rpcEXCEPT
    else:
      if verbose: print '\tfinished'
      if res is None:
        evt.result = []
      else:
        evt.result = res
      evt.rpcStatus = _rpcDONE

  # release the lock and let everyone waiting on us know that we
  #   finished
  if verbose: print 'unlock'
  evt.rpcStatusLock.unlock()
  if verbose: print 'wake'
  evt.rpcCondVar.wakeAll()
  if verbose: print '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-'
    
    
  
class rpcEvent(QCustomEvent):
  """ a QCustomEvent subclass that we'll use to actually pass
    along the event

  """
  def __init__(self,method,params):
    """ constructor

    **Arguments**

      - method: name of the method to be invoked

      - params: tuple of parameters to be passed to the method

    """
    QCustomEvent.__init__(self,_rpcEventID)
    self._method = method
    self._params = params
    self.rpcStatus = _rpcPENDING
    self.rpcStatusLock = QMutex()
    self.rpcCondVar = QWaitCondition()
    self.result = ""
  def method(self):
    """ getter for the method attr """
    return self._method
  def params(self):
    """ getter for the params attr """
    return self._params

class rpcClient:
  """ client class that receives xmlrpc events

  """
  def __init__(self,app,win,awaitResponse=1,verbose=1):
    self._verbose = verbose
    if self._verbose: print 'rpc client'
    self._app = app
    self._win = win
    self._stop = 0
    self._responsive=awaitResponse
  def _dispatch(self,method,params):
    """ dispatcher method, the _SimpleXMLRPCServer_ calls this every
        time a call is made

    """
    if method == 'quit':
      self._stop = 1
    if self._verbose: print 'client: %s %s'%(method,params)
    evt = rpcEvent(method,params)
    evt.rpcStatusLock.lock()
    evt.rpcStatus = _rpcPENDING
    evt.rpcStatusLock.unlock()

    QThread.postEvent(self._win,evt)
    if self._responsive:
      evt.rpcCondVar.wait()
    evt.rpcStatusLock.lock()
    if evt.rpcStatus == _rpcEXCEPT:
      # The GUI threw an exception, release the status lock
      #  and re-raise the exception
      evt.rpcStatusLock.unlock()
      raise RPCError,'%s: %s'%(str(evt.result[0]),str(evt.result[1]))
    else:
      if self._responsive:
        s = evt.result
      else:
        s = 'ok'
      evt.rpcStatusLock.unlock()
      return s


class ServThread(QThread):
  """ a thread which runs the xmlrpc server

  """
  def  __init__(self,instance,port=8234,nToTry=10,server=None,verbose=2):
    """ constructor

    **Arguments**

      - instance: an _rpcClient_ instance (or something equivalent) to
        process rpc calls

      - port: (optional) initial port for the server

      - nToTry: (optional) maximum number of ports that will be tried
        as we try to start the server

      - server: (optional) _SimpleXMLRPCServer_ (or something
        equivalent) instance.

      - verbose: (optional) sets the verbosity level
      

    """
    QThread.__init__(self)
    self.instance = instance
    self._verbose = verbose

    import SimpleXMLRPCServer
    i = 0
    while i < nToTry and server is None:
      try:
        server = SimpleXMLRPCServer.SimpleXMLRPCServer(("localhost",port+i),
                                                       logRequests=0)
      except:
        i+=1
        server = None
    if server is None:  raise NoServer,'could not open RPC server'
    if verbose>0: print 'xmlrpc server active on port %d'%(port+i)
    self.serverPort = port+i
    self.server = server
    self.server.register_instance(instance)

  def run(self):
    """ launches the thread """
    while( not self.instance._stop ):
      if self._verbose: print 'running'
      self.server.handle_request()
    if self._verbose: print 'done'

if __name__ == '__main__':
  import sys
  class handler:
    def __init__(self):
      self._stop = 0
    def _dispatch(self,method,params):
      print 'dispatch:',method,params
      if method == 'quit':
        self._stop = 1
      return 'ok: %s'%method

  h = handler()
  t = ServThread(h)
  print 'start serv thread'
  t.start()
  doGui = 1
  if doGui:
    a = QApplication(sys.argv)
    w = QMainWindow()
    a.setMainWidget(w)
    w.show()
    print 'start exec loop'
    a.exec_loop()
  print 'wait'
  t.wait()
  #t2.wait()
  
