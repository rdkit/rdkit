# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" Mixin GUI launcher code """
import qt
import os.path
from rdkit import qtGui
from rdkit.qtGui import qtUtils

from rdkit import RDConfig
lookHere = [os.path.join(RDConfig.RDCodeDir,'qtGui/VisProg'),
            os.path.join(RDConfig.RDCodeDir,'qtGui/MixIns'),
            '.']

def StandardGui(args,**kwargs):
  """ launches a standard mixin GUI

  """
  from rdkit.mixins import Loader
  loadOrder = []
  loadedMixIns=Loader.LoadMixIns(args,loadOrder,
                                 lookHere)
  window = qtGui.GuiBase.GuiBase(name=args[-1],**kwargs)
  Loader.InitMixIns(window,loadedMixIns)
  window.finalizeMenus()
  return window

def Launcher(klass,*args,**kwargs):
  app = qt.QApplication([klass.__name__])
  widget = klass(*args,**kwargs)
  widget.show()
  app.setMainWidget(widget)
  return app,widget


def Usage():
  """ displays a usage message and exits

  """
  import sys
  msg = """
Usage: Gui [-s] <list of mixin names>

Arguments:
  -s: start an xmlrpc server

  """
  print msg
  sys.exit(-1)

if __name__ == '__main__':
  import sys,getopt
  try:
    args,extras = getopt.getopt(sys.argv[1:],'us')
  except:
    Usage()
  doUnits = 0
  startServer = 0
  for arg,val in args:
    if arg == '-u':
      doUnits = 1
    elif arg == '-s':
      startServer = 1
  if len(extras) > 0:
    args = extras
  else:
    args = ['CompositeBuilder','MolCluster','MolDescriptors','ClusterInteract','MolSimilarity']
  app = qt.QApplication(sys.argv)
  qt.QApplication.setOverrideCursor(qt.QCursor(qt.Qt.WaitCursor))
  splash,pix = qtUtils.ShowSplashScreen(app)
  if not doUnits:
    win = StandardGui(args,startServer=startServer)
  else:
    win = UnitTestGui(args)
  splash.hide()
  win.show()
  win.updateGeometry()
  win.setCaption('RDKit Gui')
  qt.QApplication.restoreOverrideCursor()
  win.setActiveWindow()
  win.raiseW()
  app.setMainWidget(win)
  app.exec_loop()

