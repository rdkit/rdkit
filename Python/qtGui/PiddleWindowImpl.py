# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" implementation bits for PiddleWindows

 this class handles interaction

"""    
from pyRDKit import RDConfig
import qt
import qtcanvas

from pyRDKit.qtGui.forms.PiddleWindow import PiddleWindow as _Form
from pyRDKit.qtGui import qtUtils
from pyRDKit.sping.Qt.pidQt import QtCanvas as Canvas
import os.path
from pyRDKit.Logger import Logger

class PiddleCanvasView(qtcanvas.QCanvasView):
  """  The actual canvas view which is displayed to the user

    **Notes**
      - we maintain a list of event handlers for left and right
        clicks, each time one of these events comes in, those handlers
        are called in the order in which they were inserted

  """
  def __init__(self,*args,**kwargs):
    qtcanvas.QCanvasView.__init__(self,*args)
    self._leftClickHandlers = []
    self._rightClickHandlers = []
    self._leftDblClickHandlers = []
    self._rightDblClickHandlers = []
    self._tformMat = qt.QWMatrix()
    self._tformMat.setMatrix(1,0.0,0.0,1,0.0,0.0)
    self.setWorldMatrix(self._tformMat)
    self.rescaleOnResize=True
    self._sz=None

  def addLeftHandler(self,fn):
    """ adds a left click handler.

    **Arguments**:

      - fn: a function which takes a two arguments:

        1) view: this object

        2) evt: the QMouseEvent which came in
        

    """
    self._leftClickHandlers.append(fn)
  def addRightHandler(self,fn):
    """ adds a right click handler.

    **Arguments**:

      - fn: a function which takes a two arguments:

        1) view: this object

        2) evt: the QMouseEvent which came in
        

    """
    self._rightClickHandlers.append(fn)

  def addLeftDblHandler(self,fn):
    """ adds a left double click handler.

    **Arguments**:

      - fn: a function which takes a two arguments:

        1) view: this object

        2) evt: the QMouseEvent which came in
        

    """
    self._leftDblClickHandlers.append(fn)

  def addRightDblHandler(self,fn):
    """ adds a right double click handler.

    **Arguments**:

      - fn: a function which takes a two arguments:

        1) view: this object

        2) evt: the QMouseEvent which came in
        

    """
    self._rightDblClickHandlers.append(fn)
    
  def __leftMouseEvent(self,evt,dbl=0):
    """ INTERNAL USE ONLY

        process left mouse events by calling each of our registered
        handlers

    **Arguments**

      - evt: the mouse event which came in
      
    """
    if not dbl:
      fns = self._leftClickHandlers
    else:
      fns = self._leftDblClickHandlers
    for fn in fns:
      fn(self,evt)
    
  def __rightMouseEvent(self,evt,dbl=0):
    """ INTERNAL USE ONLY

        process left mouse events by calling each of our registered
        handlers

    **Arguments**

      - evt: the mouse event which came in
      
    """
    if not dbl:
      fns = self._rightClickHandlers
    else:
      fns = self._rightDblClickHandlers
    for fn in fns:
      fn(self,evt)

    
  def contentsMousePressEvent(self,evt):
    """ callback for when a mouse press happens in our contents

     The event is dispatched to the appropriate handler and then the
     canvas is updated.

    **Arguments**

      - evt: the mouse event which came in

    """
    but = evt.button()
    if but == qt.Qt.LeftButton:
      self.__leftMouseEvent(evt)
    elif but == qt.Qt.RightButton:
      self.__rightMouseEvent(evt)
    self.canvas().update()

  def contentsMouseDoubleClickEvent(self,evt):
    """ callback for when a mouse press happens in our contents

     The event is dispatched to the appropriate handler and then the
     canvas is updated.

    **Arguments**

      - evt: the mouse event which came in

    """
    but = evt.button()
    if but == qt.Qt.LeftButton:
      self.__leftMouseEvent(evt,dbl=1)
    elif but == qt.Qt.RightButton:
      self.__rightMouseEvent(evt,dbl=1)
    self.canvas().update()

  def resizeEvent(self,evt):
    # delegate, then resize the canvas
    if self._sz is None:
      self._sz = evt.size()
      return
    qtcanvas.QCanvasView.resizeEvent(self,evt)
    if self.rescaleOnResize:
      sz = evt.size()
      oldSz = evt.oldSize()
      xs = float(sz.width())/oldSz.width()
      ys = float(sz.height())/oldSz.height()
      self._tformMat.scale(xs,ys)
      self.setWorldMatrix(self._tformMat)
    canv = self.canvas()
    if canv:
      canv.update()

  def reset(self):
    self._tformMat.reset()
    self.setWorldMatrix(self._tformMat)

    


class PiddleWindow(_Form):
  """ the window which is used to contain PiddleCanvasViews

   **Important Attributes**

     - canvas: the canvas into which we draw (actually a
       _Logger.Logger_ instance around a QtCanvas

     - view: our PiddleCanvasView

     - qCanvas: INTERNAL USE ONLY the actual QCanvas associated with
       the QtCanvas (*oh what tangled webs we weave*)

  """
  def __init__(self,*args,**kwargs):
    _Form.__init__(self,*args,**kwargs)
    # set up the canvas view
    self.setCentralWidget(qt.QWidget(self,"qt_central_widget"))
    self._mainLayout = qt.QVBoxLayout(self.centralWidget(),2,2)
    self._qCanvasView = PiddleCanvasView(self.centralWidget())
    self._mainLayout.addWidget(self._qCanvasView)

    self.statusBar().setSizeGripEnabled(0)
    id = self.fileMenu.insertItem(self.trUtf8('&Export'),-1,0)
    self.fileMenu.connectItem(id,self.exportCanvas)
    self._dir = '.'
  
  def canvas(self):
    """ getter for our _canvas_ attribute """
    return self._canvas

  def view(self):
    """ getter for our _view_ attribute """
    return self._qCanvasView

  def resizeCanvas(self,canvSize):
    """ resizes our canvas

    **Arguments**

      - canvSize: a 2-sequence with the canvas dimensions

    **Notes**

      - this actually destroys the current canvas and creates a new one.

    """
    self._qCanvas = qtcanvas.QCanvas(canvSize[0],canvSize[1])
    self._qCanvasView.setCanvas(self._qCanvas)
    self._canvas = Logger.Logger(Canvas,destCanvas=self._qCanvas,size=canvSize,
                                 loggerFlushCommand='clear',loggerIgnore=['size'])

  def exportCanvas(self,destCanvas=None):
    """ exports our canvas using the _Logger_ features

     **Arguments**

       - destCanvas: (optional) The Piddle/Sping Canvas onto which the
         _Logger_ should replay.  If this is not provided, we'll get
         one from _self.setupPiddleExportCanvas_
         
    """
    canvas = self.setupPiddleExportCanvas()
    if canvas:
      Logger.replay(self._canvas,canvas)
      canvas.save()
    
  def setupPiddleExportCanvas(self,size=None,extension=''):
    """ sets up and returns a Piddle/Sping canvas

    **Arguments**

      - size: (optional) the size of the target canvas.  If this is
        not provided, we'll duplicate the size of our existing canvas.

      - extension: (optional) the type of extension to allow for
        output canvases.  If this is not provided, we'll construct a
        list of all available extensions.

    **Returns**

      a fresh Piddle/Sping canvas of the appropriate type

    """
    if size is None:
      size = self._canvas.getSize()
    resCanv,fileName = qtUtils.GetPiddleCanvas(size,extension=extension,dir=self._dir)
    # FIX: need some error checking here
    if resCanv and fileName:
      self._dir = os.path.dirname(fileName)
    return resCanv
  #
  # Slots
  #
  def fileClose(self):
    """ callback for File->Close

    """
    self.hide()

  def resizeEvent(self,evt):
    pass
if __name__ == '__main__':
  from pyRDKit.qtGui import Gui

  app,widg = Gui.Launcher(PiddleWindow,None)
  widg.resize(640,480)
  widg.resizeCanvas((640,480))
  app.exec_loop()
  widg.destroy(1)
