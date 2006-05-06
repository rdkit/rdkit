# $Id$
#
#  Copyright (C) 2002-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
import RDConfig
from qt import *
import os.path
from sping.utils import availableCanvases
import sys

logoImageData=[
"16 16 6 1",
". c None",
"# c #1414ff",
"a c #143cff",
"b c #6464ff",
"c c #8c8cff",
"d c #b4b4ff",
"................",
".....######.....",
"....#aaaaaa#....",
"...#abbbbbba#...",
"..#abbbbbbbba#..",
".#abbbccccbbba#.",
".#abbccccccbba#.",
".#abbccddccbba#.",
".#abbccddccbba#.",
".#abbccccccbba#.",
".#abbbccccbbba#.",
"..#abbbbbbbba#..",
"...#abbbbbba#...",
"....#aaaaaa#....",
".....######.....",
"................"
]


# used to assign a unique ID to every item on the canvas
_objsSoFar=0

logger = None
if 0:
  import logging
  def _initLogger():
    global logger
    logger = logging.getLogger('qtUtils')
    hdlr = logging.StreamHandler(sys.stderr)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
else:
  import RDLogger as logging
  def _initLogger():
    global logger
    logger = logging.logger()
_initLogger()

def error(msg,parent=None,caption="Error",**kwargs):
  """ displays an error message """
  logger.error(msg,**kwargs)
  mb = QMessageBox.critical(parent,caption,"ERROR: %s"%(msg),
                            QMessageBox.Ok,QMessageBox.NoButton,QMessageBox.NoButton)
  return mb

def warning(msg,parent=None,caption="Warning",**kwargs):
  """ displays a warning message """
  logger.warning(msg,**kwargs)
  mb = QMessageBox.warning(parent,caption,"Warning: %s"%(msg),
                           QMessageBox.Ok,QMessageBox.NoButton,QMessageBox.NoButton)
  return mb
def information(msg,parent=None,caption="Note",**kwargs):
  """ displays an informational message """
  logger.info(msg,**kwargs)
  mb = QMessageBox.information(parent,caption,"%s"%(msg),
                           QMessageBox.Ok,QMessageBox.NoButton,QMessageBox.NoButton)
  return mb

def infoPrompt(msg,parent=None,caption="Note",**kwargs):
  """ displays a warning message """
  logger.info(msg,**kwargs)
  mb = QMessageBox.information(parent,caption,"%s"%(msg),
                           QMessageBox.Yes,QMessageBox.No,QMessageBox.NoButton)
  return mb

def GetSplashPixmap(fileN=None):
  if fileN is None:
    fileN = os.path.join(RDConfig.RDImageDir,"RD-Splash.jpg")
  try:
    open(fileN,'r')
    pix = QPixmap(fileN)
  except:
    import traceback
    traceback.print_exc()
    pix = None
  return pix

def ShowSplashScreen(app,fileN=None):
  splash = QLabel(None,"RDGui Splash",
                  Qt.WStyle_Customize|Qt.WStyle_NoBorder|Qt.WStyle_Tool)
  pix = GetSplashPixmap(fileN=fileN)
  if pix:
    splash.setPixmap(pix)
    d = app.desktop()
    splash.move((d.width()-pix.width())/2,(d.height()-pix.height())/2)
    splash.setCaption("RDGui Splash")
    splash.show()
    app.processEvents()

  return splash,pix
  
def GetPiddleCanvas(size,extension='',dir=''):
  """ sets up and returns a Piddle/Sping canvas

    **Arguments**

      - size: the size of the target canvas.

      - extension: (optional) the type of extension to allow for
        output canvases.  If this is not provided, we'll construct a
        list of all available extensions.

    **Returns**

      a 2-tuple:

        - a fresh Piddle/Sping canvas of the appropriate type

        - the filename  

  """
  if extension == '' or not availableCanvases.has_key(extension):
    typeList = []
    for key in availableCanvases.keys():
      tStr = '%s files (*.%s)'%(availableCanvases[key][2],key)
      typeList.append(tStr)
    typeList = ';;'.join(typeList)
    fileName = str(QFileDialog.getSaveFileName(dir,typeList))
    extension = fileName.split('.')[-1]
  else:
    moduleName,canvasName,typeName = availableCanvases[extension]
    prompt = 'Choose %s FileName'%typeName
    typeList = '%s Files (*.%s);;All Files (*.*)'%(typeName,extension)
    fileName = str(QFileDialog.getSaveFileName(dir,typeList))
  if not fileName:
    return None,''
  # FIX: need some error checking here
  moduleName,canvasName,typeName = availableCanvases[extension]
  exec('from sping.%s import %s'%(moduleName,canvasName))
  exec('resCanv = %s(size,"%s")'%(canvasName,fileName))
  return resCanv,fileName

class ProgressStop(StopIteration):
  pass

class ProgressDialogHolder:
  def __init__(self,label,endPt):
    self.dlg = QProgressDialog(label,'Cancel',endPt)
    self.dlg.setLabelText(label)
    self.dlg.setCaption(label)
    image0 = QPixmap(logoImageData)
    self.dlg.setIcon(image0)
  def __call__(self,val):
    self.dlg.setProgress(val)
    qApp.processEvents()
    if self.dlg.wasCancelled():
      raise ProgressStop,'aborted'
    
def safeDestroyWidget(widget):
  widget.reparent(QWidget(),0,QPoint())

def getSaveFilenameAndFilter(dir,filter,mode=QFileDialog.AnyFile):
  dlg = QFileDialog(dir,filter)
  dlg.setMode(mode)
  r = dlg.exec_loop()
  fName = str(dlg.selectedFile())
  filter = str(dlg.selectedFilter())
  return fName,filter
  
