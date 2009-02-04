# $Id$
#
# Copyright (C) 2005,2006 Rational Discovery LLC
#  All Rights Reserved
#
_version='0.9.47'

import os,sys
from utils import comhack
if hasattr(sys,'frozen') and sys.frozen==1:
  dirPath=os.path.join(os.path.dirname(sys.executable),'support','gen_py')
  try:
    comhack.set_gen_path(dirPath)
  except OSError:
    print >>sys.stderr,'Failed using default directory for gen_py files.\nUsing a temporary directory instead.'
    dirPath=os.path.join(os.environ.get('TEMP','.'),'support','gen_py')
    try:
      comhack.set_gen_path(dirPath)
    except:
      import traceback
      traceback.print_exc()
      print >>sys.stderr,'Cannot open temp directory for gen_py files.'
      sys.exit(-1)
  os.environ['RDBASE']=os.path.dirname(sys.executable)
  os.environ['PYTHONPATH']=''

from qtGui.Search3D import LocalConfig
if hasattr(sys,'frozen') and sys.frozen==1:
  LocalConfig.applicationName='RDPharm3D'
  LocalConfig.neighborhoodRadius=7.5
  LocalConfig.logFilename='RDPharm3D_log.txt'
from qt import *
# we put this up here so that the splash can display while everything imports

_splashMessage="""
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  %s version %s

  Copyright (C) 2005-2007 Greg Landrum

  This open-source software is part of RDKit (www.rdkit.org).
  Please see the file license.txt for information about the
  license.

  This software is copyrighted.  The software may not be copied,
  reproduced, translated or reduced to any electronic medium or
  machine-readable form without the prior written consent of
  Rational Discovery LLC.
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
"""%(LocalConfig.applicationName,_version)


import os,sys
if not os.environ.has_key('RDBASE') or not os.environ['RDBASE']:
  os.environ['RDBASE']=os.getcwd()
import RDConfig
#from utils import chemdraw
splashFilename = os.path.join(RDConfig.RDDocsDir,'Programs',LocalConfig.applicationName,
                              '%s-Splash.jpg'%LocalConfig.applicationName)
if __name__=='__main__':
  import sys,time
  from qtGui import Gui,GuiBase,qtUtils
  app = QApplication(sys.argv)
  QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
  print _splashMessage
  splash,pix = qtUtils.ShowSplashScreen(app,fileN=splashFilename)
  t1 = time.time()

from qtGui import GuiBase,qtUtils

enableDbWidget=False

import Chem
from qtGui.Search3D.Search3DWidget import Search3DWidget
from qtGui.Search3D import LocalConfig
from PIL import Image,ImageDraw,ImageFont
# these are here to help the Installer work:
import PIL.ImageFile
import PIL.JpegImagePlugin
import PIL.BmpImagePlugin
import PIL.PngImagePlugin


if hasattr(sys,'frozen') and sys.frozen==1:
  if not os.environ.has_key('RDPHARM3D_EXPERT'):
    LocalConfig.allowExpertMode=False
  else:
    LocalConfig.allowExpertMode=True

#logPath=os.environ.get('TEMP','.')
logPath=os.getcwd()
LocalConfig.logFilename = os.path.join(logPath,LocalConfig.logFilename)
qtUtils.logger.setLevel(qtUtils.logging.DEBUG)
qtUtils.logging.AttachFileToLog('rdApp.error',LocalConfig.logFilename)
qtUtils.logging.AttachFileToLog('rdApp.warning',LocalConfig.logFilename)
qtUtils.logging.AttachFileToLog('rdApp.debug',LocalConfig.logFilename)
qtUtils.logging.AttachFileToLog('rdApp.info',LocalConfig.logFilename)
class Search3DApp(GuiBase.GuiBase):
  def __init__(self,*args,**kwargs):
    kwargs['name']  = LocalConfig.applicationName
    kwargs['includeLogo']  = False
    GuiBase.GuiBase.__init__(self,*args,**kwargs)
    self._initHelp()
    self._dir = '.'
    self._filename = ''
    self.setCentralWidget(QWidget(self,'central_widget'))
    self.enableDbWidget=enableDbWidget

    if not self.centralWidget():
      self.setCentralWidget(QWidget(self,'central_widget'))
      
    self.searchWidget = Search3DWidget(self.centralWidget())
    layout = QHBoxLayout(self.centralWidget(),1,1,'centrallayout')
    layout.addWidget(self.searchWidget)
    self.statusBar().setSizeGripEnabled(True)

    self._actionMenu = QPopupMenu(self)
    self.menuBar().insertItem(self.trUtf8('&Actions'),self._actionMenu)
    self.addMenuItems()
    self.finalizeMenus()

    self.enableExpertMode(state=self.expertModeEnabled())
    
  def addMenuItems(self):
    self._mbExportMenu = QPopupMenu(self)
    self._fileOpenId=self._fileMenu.insertItem(self.trUtf8('&Open'),self.openState)
    self._fileSaveId=self._fileMenu.insertItem(self.trUtf8('&Save'),self.saveState)
    self._fileSaveAsId=self._fileMenu.insertItem(self.trUtf8('Save &As'),self.saveStateAs)
    self._fileExportId=self._fileMenu.insertItem(self.trUtf8('E&xport'),self._mbExportMenu,0)
    self._mbExportMenu.insertItem(self.trUtf8('&All Evaluated Molecules'),self.exportAllHits)
    self._mbExportMenu.insertItem(self.trUtf8('&Checked &Alignments'),self.exportEmbeds)

    self._fileLoadParamsId=self._fileMenu.insertItem(self.trUtf8('&Load Params'),
                                                     self.loadParams)

    self._actionMenu.insertItem(self.trUtf8('&Cleanup Viewer'),
                                self.searchWidget.cleanupViewer)
    self._actionMenu.insertItem(self.trUtf8('&Show Hydrogen Bonds'),
                                self.highlightHBonds)
    self._actionMenu.insertItem(self.trUtf8('&Load SD File'),
                                self.searchWidget.loadSdFileAction)
    self._actionMenu.insertItem(self.trUtf8('&Grab From ChemDraw'),
                                self.searchWidget.cdxGrab)

    self._mbAdvancedActionMenu = QPopupMenu(self)
    self._actionMenu.insertItem(self.trUtf8('&Advanced'),self._mbAdvancedActionMenu)
    id=self._mbAdvancedActionMenu.insertItem(self.trUtf8('&Expert Mode'),
                                             self.enableExpertMode)
    self._expertModeMenuId=id
    state = LocalConfig.allowExpertMode
    self._mbAdvancedActionMenu.setItemChecked(self._expertModeMenuId,state)
    if not state:
      self._mbAdvancedActionMenu.setDisabled(True)

    id=self._mbAdvancedActionMenu.insertItem(self.trUtf8('&Add Hs?'),
                                             self.enableAddHs)
    self._enableAddHsMenuId=id
    self._mbAdvancedActionMenu.setItemChecked(self._enableAddHsMenuId,False)

  def expertModeEnabled(self):
    if hasattr(self,'_mbAdvancedActionMenu'):
      return self._mbAdvancedActionMenu.isItemChecked(self._expertModeMenuId)
    else:
      return False
  def enableExpertMode(self,id=-1,state=None):
    if state is None:
      state = not self.expertModeEnabled()
    self._mbAdvancedActionMenu.setItemChecked(self._expertModeMenuId,state)
    if state and hasattr(self,'_mbAdvancedActionMenu'):
      self._mbAdvancedActionMenu.setEnabled(True)
    self.searchWidget.enableExpertMode(state)
    self._fileMenu.setItemEnabled(self._fileLoadParamsId,state)
    
  def enableAddHs(self,id=-1,state=None):
    if state is None:
      state = not self._mbAdvancedActionMenu.isItemChecked(self._enableAddHsMenuId)
    self._mbAdvancedActionMenu.setItemChecked(self._enableAddHsMenuId,state)
    self.searchWidget.addHs=state
    
  def _initHelp(self):
    from StringIO import StringIO
    self._aboutWin = QDialog(None)
    self._aboutWin.setCaption("About %s"%LocalConfig.applicationName)
    vout = QVBoxLayout(self._aboutWin)
    lab = QLabel(self._aboutWin,LocalConfig.applicationName,
                 Qt.WStyle_Customize|Qt.WStyle_NoBorder|Qt.WStyle_Tool)

    img = Image.open(splashFilename)
    d = ImageDraw.Draw(img)
    try:
      fnt = ImageFont.truetype('arial.ttf',18)
    except IOError:
      fnt = None
    if fnt:
      sz = img.size
      pos = sz[0]-120,sz[1]-25
      d.text(pos,'Version: %s'%_version,font=fnt,fill=(0,0,0))
      del d
      sio = StringIO()
      img.save(sio,format='bmp')
      pix = QPixmap()
      pix.loadFromData(sio.getvalue())
      lab.setPixmap(pix)
      vout.addWidget(lab)

    hout = QHBoxLayout()
    vout.addLayout(hout)
    butt = QPushButton(self._aboutWin,"Ok")
    hout.addWidget(butt)
    butt.setText("OK")
    self.connect(butt,SIGNAL('clicked()'),self._aboutWin.accept)
    self._aboutWin.adjustSize()
    
  def setDefaultDir(self,dir=None):
    if dir is None:
      dir = os.getcwd()
      
    self._dir=dir

  def finalizeMenus(self):
    """ finalizes our menu bar

    """
    GuiBase.GuiBase.finalizeMenus(self)
    self.menuBar().removeItem(self._viewMenu._id)
    self._viewMenu = None
    if hasattr(sys,'frozen') and sys.frozen==1:
      self._fileMenu.removeItem(self._fileMenu._pyShellId)
      self.menuBar().removeItem(self._editMenu._id)
    id=self._helpMenu.insertItem(self.trUtf8('&Help'),
                                 self.launchHelp,Qt.CTRL+Qt.Key_F1,-1,0)
    if os.environ.get('RDPHARM3D_HOWTOLOCATION',''):
      self._helpMenu.insertItem(self.trUtf8('HowTos and &FAQs'),
                                self.launchHowTo,0,-1,1)

  def launchHelp(self):
    if os.environ.has_key('RDPHARM3D_HELPLOCATION'):
      os.startfile(os.environ.get('RDPHARM3D_HELPLOCATION'))
    else:
      os.startfile(os.path.join(RDConfig.RDDocsDir,"Programs",LocalConfig.applicationName,
                                "%s.chm"%LocalConfig.applicationName))
  def launchHowTo(self):
    if os.environ.has_key('RDPHARM3D_HOWTOLOCATION'):
      os.startfile(os.environ.get('RDPHARM3D_HOWTOLOCATION'))

  def fileClose(self):
    """ callback for File->Close

    """
    self.hide()


  def aboutBox(self):
    self._aboutWin.exec_loop()



  def exportFile(self,mols,filename='',filter=''):
    if not mols:
      return
    if not filename:
      filename,filter = qtUtils.getSaveFilenameAndFilter(self._dir,
                                                       'SD Files (*.sdf);;SMILES Tables (*.txt);;CSV Files (*.csv);;TDT Files (*.tdt)')
    if not filename:
      return

    splitN=filename.split('.')
    if len(splitN)>1:
      ext = splitN[-1].lower()
    else:
      ext = ''
    if ext not in ('sdf','tdt','txt','csv') and filter:
      filter = filter.split(' ')[0].lower()
      if filter=='sd':
        ext='sdf'
      elif filter=='smiles':
        ext='txt'
      elif filter=='csv':
        ext='csv'
      elif filter=='tdt':
        ext='tdt'
      if ext:  
        filename += '.'+ext

    if os.path.exists(filename):
      name = os.path.basename(filename)
      confirm = QMessageBox.warning(self,
                              self.trUtf8("Overwrite file?"),
                              self.trUtf8("File %s already exists, should it be overwritten?"%(name)),
                              1,2)
      if confirm != 1:
        return

    if ext=='sdf':
      writer = Chem.SDWriter(filename)
    elif ext=='tdt':
      writer = Chem.TDTWriter(filename)
      writer.SetWriteNames(True)
    elif ext=='txt':
      writer = Chem.SmilesWriter(filename,delimiter='\t',includeHeader=True)
    elif ext=='csv':
      writer = Chem.SmilesWriter(filename,delimiter=',',includeHeader=True)
    else:
      qtUtils.error('unrecognized extension:\n%s'%ext)
      return
    self.setDefaultDir(os.path.dirname(filename))

    for m in mols:
      qtUtils.logger.debug(m.GetProp("_Name"))
      writer.write(m)
      writer.flush()

  def exportAllHits(self,id,filename=''):
    if not hasattr(self.searchWidget,'molAlignList'):
      return
    hits = self.searchWidget.molAlignList.getMols(attr='getMol2d')
    self.exportFile(hits,filename=filename)
    
  def exportEmbeds(self,id,filename=''):
    if not hasattr(self.searchWidget,'molAlignList'):
      return
    hits = self.searchWidget.molAlignList.getCheckedEmbeds()
    self.exportFile(hits,filename=filename)

  def highlightHBonds(self,id):
    widg = self.searchWidget
    if not widg.proteinName or \
       not hasattr(widg,'mol') or not widg.mol:
      return
    widg.displayHBonds()

  def saveState(self,id=-1,filename=None):
    if filename is None:
      filename = self._filename
    if not filename:
      self.saveStateAs()
      return

    outData = self.searchWidget.toString()
    try:
      outF = file(filename,'wb+')
    except:
      qtUtils.error('could not open file %s for writing'%filename,exc_info=True)
    else:
      outF.write(outData)
      outF.close()
      self._filename = filename

  def saveStateAs(self,id=-1):
    filename = str(QFileDialog.getSaveFileName(self._dir,
                                               '%s files (*.p3d);;All files (*.*)'%LocalConfig.applicationName))
    if not filename:
      return
    else:
      self.setDefaultDir(os.path.dirname(filename))
      if filename.find('.p3d')<0:
        filename = filename + '.p3d'

      if os.path.exists(filename):
        name = os.path.basename(filename)
        confirm = QMessageBox.warning(self,
                                self.trUtf8("Overwrite file?"),
                                self.trUtf8("File %s already exists, should it be overwritten?"%(name)),
                              1,2)
        if confirm != 1:
          return


      self.saveState(filename=filename)
                     
  def openState(self,id=-1,filename=None):
    if filename is None:
      filename = str(QFileDialog.getOpenFileName(self._dir,
                                                 '%s files (*.p3d);;All files (*.*)'%LocalConfig.applicationName))
      if filename:
        self.setDefaultDir(os.path.dirname(filename))
    if not filename:
      return False
    else:
      try:
        inF = open(filename,'rb')
      except:
        qtUtils.error('could not open file %s for reading'%filename,exc_info=True)
        return False
      else:
        inD = inF.read()
        inF.close()
        self.searchWidget.fromString(inD)
        #self._filename = filename

        self._fileMenu.setItemEnabled(self._fileSaveId,True)
        self._fileMenu.setItemEnabled(self._fileSaveAsId,True)
        self._fileMenu.setItemEnabled(self._fileExportId,True)
        self.searchWidget.grabAtomsButton.setEnabled(True)

        return True

  def loadParams(self,id=-1,filename=None):
    if filename is None:
      filename = str(QFileDialog.getOpenFileName(self._dir,
                                                 'Python files (*.py);;All files (*.*)'))
      if filename:
        self.setDefaultDir(os.path.dirname(filename))
    if not filename:
      return False
    else:
      try:
        open(filename,'r')
      except:
        qtUtils.error('could not open file %s for reading'%filename,exc_info=True)
        return False
      else:
        try:
          execfile(filename,LocalConfig.__dict__)
        except:
          qtUtils.error('problems importing parameter file %s'%filename,exc_info=True)
          return False
        return True

    
if __name__=='__main__':
  Search3DApp.fileClose = Search3DApp.closeEvent
  window = Search3DApp()

  t2 = time.time()
  # make sure the splash stays visible for a minimum
  # amount of time:
  timeRemain = 1 - (t2 - t1)
  if timeRemain>0:
    time.sleep(timeRemain)
  splash.hide()

  window.show()
  window.updateGeometry()
  window.adjustSize()
  window.setActiveWindow()
  window.raiseW()

  window.setDefaultDir()
  if len(sys.argv)>1:
    window.openState(filename=sys.argv[1])
    window.enableExpertMode(state=False)
    QApplication.restoreOverrideCursor()
  elif hasattr(sys,'frozen') and sys.frozen==1:
    QApplication.restoreOverrideCursor()
    if not LocalConfig.allowExpertMode:
      window.enableExpertMode(state=False)
      ok = window.openState()
      if not ok:
        window._fileMenu.setItemEnabled(window._fileSaveId,False)
        window._fileMenu.setItemEnabled(window._fileSaveAsId,False)
        window._fileMenu.setItemEnabled(window._fileExportId,False)
        window.searchWidget.grabAtomsButton.setEnabled(False)
  else:
    QApplication.restoreOverrideCursor()


  # Issue314:
  # This is specific to Sunesis, who want different
  # tab names for their version:
  if hasattr(sys,'frozen') and sys.frozen==1:
    window.searchWidget.tabWidget.setTabLabel(window.searchWidget.pcophorePage,
                                              'Alignment Criteria')
    window.searchWidget.tabWidget.setTabLabel(window.searchWidget.detailsPage,
                                              'Advanced Search')
    window.searchWidget.tabWidget.setTabLabel(window.searchWidget.resultsPage,
                                              'Molecule Evaluation')
    
  app.exec_loop()
