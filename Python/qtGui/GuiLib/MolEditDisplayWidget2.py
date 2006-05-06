# $Id$
#
#  Copyright (C) 2003-2005  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" implementation bits for MolEditDisplayWidgets

"""    
import RDConfig
from qt import *
from qtGui.GuiLib.forms.MolEditDisplayWidget2 import MolEditDisplayWidget as _Form
from qtGui import qtUtils
import cStringIO as StringIO
from PIL import Image
from utils import PilTools
import os.path,types,re
try:
  from utils import chemdraw
except ImportError:
  chemdraw = None
else:
  try:
    chemdraw.StartChemDraw
  except AttributeError:
    chemdraw = None
from qtGui.GuiLib import MolCanvas

fileFormats={
  'mol':('Mol Files (*.mol)','chemical/x-mdl-molfile'),
  'sdf':('SD Files (*.sdf)','chemical/x-mdl-sdfile'),
  'mol2':('Mol2 Files (*.mol2)','chemical/msi-mol2file'),
  }

class MolEditDisplayWidget(_Form):
  """ DOC
  

  """
  def __init__(self,parent=None,initDir='',startCDX=1,allowArrows=0,logging=0):
    _Form.__init__(self,parent)
    self._dir = initDir
    self.activeMol = None
    self.molLog = []
    self.activeIdx = -1
    self.molImage = None
    self.fmt='chemical/daylight-smiles'
    self.grabFmt='chemical/daylight-smiles'

    self._molLoadFormats = []
    self._supportedDropFormats = []
    for key,val in fileFormats.iteritems():
      self._supportedDropFormats.append(key)
      self._molLoadFormats.append(val[0])
    self._molUpdateCallbacks = []
    self._nextMolCallbacks = []
    self._prevMolCallbacks = []
    self._internalUpdate=0
    self._smilesLogged=0
    self.setAcceptDrops(1)
    self._initChemdraw(startCDX)
    if allowArrows:
      self._initArrows()
    else:
      self.prevButton = None
      self.nextButton = None
      self.molNumLable = None
    self._logging = logging
    self.connect(self.smilesEdit,SIGNAL("returnPressed()"),self.smilesReturnPressed)
    self.constantUpdates = 1
    self._toSmilesConvertor = None

    molCanv=MolCanvas.MolCanvasView(self)
    molCanv.initCanvas()
    self.layout().insertWidget(0,molCanv)
    self.molCanvas=molCanv
    
  def _initArrows(self):
    layout = QHBoxLayout(None,2,4,"formatGroupLayout")
    self.prevButton = QPushButton(self.formatGroup,"prevMolButton")
    self.prevButton.setEnabled(0)

    self.prevButton.setText(self.trUtf8("Previous"))
    QToolTip.add(self.prevButton,self.trUtf8("Move backwards in the list of mol\
ecules."))
    layout.addWidget(self.prevButton)
    spacer_2 = QSpacerItem(60,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout.addItem(spacer_2)

    self.molNumLabel = QLabel(self.formatGroup,"molNumLabel")
    self.molNumLabel.setText(self.trUtf8(""))
    self.molNumLabel.setAlignment(QLabel.AlignCenter)
    QToolTip.add(self.molNumLabel,self.trUtf8("Number of the current molecule i\
n the molecule list."))
    layout.addWidget(self.molNumLabel)
    spacer_3 = QSpacerItem(60,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
    layout.addItem(spacer_3)

    self.nextButton = QPushButton(self.formatGroup,"nextButton")
    self.nextButton.setEnabled(0)
    self.nextButton.setText(self.trUtf8("Next"))
    QToolTip.add(self.nextButton,self.trUtf8("Move forwards in the list of mol\
ecules."))
    layout.addWidget(self.nextButton)
    self.formatGroup.layout().addLayout(layout)
 
    self.connect(self.prevButton,SIGNAL("clicked()"),self.previousButtonClicked)
    self.connect(self.nextButton,SIGNAL("clicked()"),self.nextButtonClicked)

    self.formatGroup.updateGeometry()
    self.updateGeometry()
    
  def _initChemdraw(self,startCDX):
    global chemdraw
    if chemdraw:
      if startCDX:
        try:
          chemdraw.StartChemDraw(visible=1)
        except:
          qtUtils.warning("Problems were encountered initializing ChemDraw",exc_info=True)
          chemdraw = None
          return
      self.GrabMol=chemdraw.CDXGrab
      self.PutMol=chemdraw.CDXDisplay
      self.fromCDXButton.setEnabled(1)
    self.insertMolUpdateCallback(self.getMolPicture)

  def getMolPicture(self,data,format):
    import Chem
    if format=='chemical/daylight-smiles':
      try:
        mol = Chem.MolFromSmiles(data)
      except:
        qtUtils.warning('Could not convert smiles %s to a molecule'%repr(data),exc_info=True)
        mol=None
    elif format=='chemical/x-mdl-molfile':
      try:
        mol = Chem.MolFromMolBlock(data)
      except:
        qtUtils.warning('Could not convert smiles %s to a molecule'%repr(data),exc_info=True)
        mol=None
    else:
      mol = None
      
    if mol:
      Chem.Kekulize(mol)
      self.molCanvas.setMol(mol)

  def getCanvas(self):
    return self.molCanvas

  def setToSmilesConvertor(self,val):
    self._toSmilesConvertor = val

  def setLogging(self,val):
    self._logging = val
  def getLogging(self):
    return self._logging

  def setFormat(self,fmt):
    self.fmt = fmt
  def getFormat(self):
    return self.fmt

  def addDropFormat(self,fmt):
    """ NOTE: this really takes a file extension, which should not include the '.'
    """
    self._supportedDropFormats.append(fmt)
  def insertMolLoadFormat(self,fmt,pos=-2):
    """ NOTE: this really takes a file extension, which should not include the '.'
    """
    if pos < -1:
      pos = len(self._molLoadFormats) + pos + 1
    self._molLoadFormats.insert(pos,fmt)
  
  def getMol(self):
    return self.activeMol
  def logMols(self,data):
    """ takes a sequence of (data,format) tuples """
    if type(data) not in (types.ListType,types.TupleType):
      data = [(data, self.getFormat())]
    self.activeIdx = self.getNumMols()
    self.molLog += list(data)
    self._updateButtonState()
  def setMol(self,data):
    self.activeMol=data
    if data and chemdraw:
      self.toCDXButton.setEnabled(1)
  def getNumMols(self):
    return len(self.molLog)
  def getActiveIdx(self):
    return self.activeIdx

  def _updateButtonState(self):
    if not self.prevButton or not self.nextButton: return
    if self.activeIdx > 0:
      self.prevButton.setEnabled(1)
    else:
      self.prevButton.setEnabled(0)
    if self.activeIdx < self.getNumMols()-1:
      self.nextButton.setEnabled(1)
    else:
      self.nextButton.setEnabled(0)
    self.molNumLabel.setText("%d / %d"%(self.activeIdx+1,self.getNumMols()))

  def nextMol(self):
    if self.smilesRadio.isChecked() and not self._smilesLogged:
      self.logMols([(self.getMol(),self.getFormat())])
      self._smilesLogged=1

    if self.activeIdx>=self.getNumMols():
      raise IndexError,'attempt to move off end of mol list'
    self.activeIdx += 1
    data,fmt = self.molLog[self.activeIdx]
    self.setMol(data)
    self.setFormat(fmt)
    if fmt=='chemical/daylight-smiles':
      self._internalUpdate=1
      self.smilesEdit.setText(data)
      self.smilesRadio.setChecked(1)
      self._internalUpdate=0
    else:
      self.smilesRadio.setChecked(0)
      
    self._updateButtonState()

  def prevMol(self):
    if self.smilesRadio.isChecked() and not self._smilesLogged:
      self.logMols([(self.getMol(),self.getFormat())])
      self._smilesLogged=1

    if self.activeIdx<1:
      raise IndexError,'attempt to move off beginning of mol list'
    self.activeIdx -= 1
    data,fmt = self.molLog[self.activeIdx]
    self.setMol(data)
    self.setFormat(fmt)
    if fmt=='chemical/daylight-smiles':
      self._internalUpdate=1
      self.smilesEdit.setText(data)
      self.smilesRadio.setChecked(1)
      self._internalUpdate=0
    else:
      self.smilesRadio.setChecked(0)
      
    self._updateButtonState()

  def getMolImage(self):
    return self.molImgData

  def clearMolUpdateCallbacks(self):
    self._molUpdateCallbacks=[]
  def insertMolUpdateCallback(self,cb,pos=-1):
    """

      the callback should take two arguments:
        1) data string
        2) data format string (a mime type)
      and should return 0 on success, anything else otherwise

    """
    if pos < -1:
      pos = len(self._molUpdateCallbacks) + pos + 1
    self._molUpdateCallbacks.insert(pos,cb)

  def clearNextMolCallbacks(self):
    self._nextMolCallbacks=[]
  def insertNextMolCallback(self,cb,pos=-1):
    """

      the callback should take one arguments:
        1) an integer with the index of the currently active mol

    """
    if pos < -1:
      pos = len(self._nextMolCallbacks) + pos + 1
    self._nextMolCallbacks.insert(pos,cb)
  def clearPrevMolCallbacks(self):
    self._prevMolCallbacks=[]
  def insertPrevMolCallback(self,cb,pos=-1):
    """

      the callback should take one arguments:
        1) an integer with the index of the currently active mol

    """
    if pos < -1:
      pos = len(self._prevMolCallbacks) + pos + 1
    self._prevMolCallbacks.insert(pos,cb)

  def _processUpdatedMol(self,data,format):
    """
     runs through the callbacks in order, stopping if one fails

     returns 0 on success, nonzero otherwise
    """
    qtUtils.logger.debug('PROCESS: %s'%str(format))
    for cb in self._molUpdateCallbacks:
      try:
        r = cb(data,format)
      except:
        r = 1
      qtUtils.logger.debug('  cb: %s %s'%(str(cb),str(r)))
      if r:
        return r
    return 0  
    
  def molUpdated(self):
    if self.smilesRadio.isChecked():
      fmt = 'chemical/daylight-smiles'
      rawD = str(self.smilesEdit.text())
      ext = 'smi'
    else:
      fileN = str(self.fileEdit.text())
      ext = fileN.split('.')[-1]
      ext.lower()
      try:
        fmt=fileFormats[ext][-1]
      except KeyError:
        qtUtils.error('Bad file format: %s'%ext)
        return 0
      rawD = open(os.path.join(self._dir,fileN),'r').read()
    if ext not in ['sdf']:
      logged=0
      ok = not self._processUpdatedMol(rawD,fmt)
      if ok:
        self.setFormat(fmt)
        self.setMol(rawD)
    elif ext=='sdf':
      logged=1
      rawD = rawD.replace('\r','')
      splitD = re.split('\$\$\$\$[\ ]*\n',rawD)
      rawD = []
      fmt = 'chemical/mdl-molfile'
      for entry in splitD:
        if entry.find('M  END')>1:
          d = entry.split('M  END')[0]+'M  END\n'
          rawD.append((d,fmt))
      if len(rawD):
        molOk = not self._processUpdatedMol(rawD[0][0],fmt)
        self.setFormat(fmt)
        self.setMol(rawD[0])
        self.logMols(rawD)
    return logged

  def grabMol(self,molFormat=None):
    res = 1
    if molFormat is None:
      molFormat = self.getFormat()
    try:
      rawD = self.GrabMol(outFormat=molFormat)
    except:
      qtUtils.warning("Problems encountered grabbing molecule",exc_info=True)
      res = 0
    else:
      if rawD:
        if molFormat == "chemical/daylight-smiles":
          self.smilesEdit.setText(rawD)
      else:
        res = 0
    return res
  def sendMol(self,molData=None,molFormat=None,clear=1):
    if molData is None:
      molData = self.getMol()
    if molFormat is None:
      molFormat = self.getFormat()
    if not molData:
      return
    try:
      self.PutMol(molData,molFormat,clear=clear)
    except:
      qtUtils.warning("Problems encountered sending molecule",exc_info=True)
      
  #-------------------
  #
  #  Slots
  # 
  #-------------------
  def fileButtonClicked(self,fileN=None):
    """ callback """
    if fileN is None:
      fileN = str(QFileDialog.getOpenFileName(self._dir,
                                              ";;".join(self._molLoadFormats)))
    if fileN:
      if self.smilesRadio.isChecked() and not self._smilesLogged:
        self.logMols([(self.getMol(),self.getFormat())])
        self._smilesLogged=1
      self.fileRadio.setChecked(1)
      self._dir,fileN = os.path.split(fileN)
      self._molName = fileN
      self.fileEdit.setText(fileN)
      logged=self.molUpdated()
      if not logged:
        self.logMols([(self.getMol(),self.getFormat())])
      
  def toCDXClick(self):
    """ callback """
    #self.molUpdated()
    self.sendMol()

  def fromCDXClick(self):
    """ callback """
    if self.smilesRadio.isChecked() and not self._smilesLogged:
      self.logMols([(self.getMol(),self.getFormat())])
      self._smilesLogged=1
    self._internalUpdate=1
    ok = self.grabMol(molFormat=self.grabFmt)
    # we go from chemdraw -> smiles
    if ok and self.grabFmt=="chemical/daylight-smiles":
      self.smilesRadio.setChecked(1)
    else:
      self.setFormat(self.grabFmt)
    self.molUpdated()
    self.toCDXButton.setEnabled(1)
    self.logMols([(self.getMol(),self.getFormat())])
    self._internalUpdate=0

  def smilesEdited(self):
    """ callback """
    if self._internalUpdate:
      # we updated the smiles internally... do not do anything
      return
  
    smi = str(self.smilesEdit.text())
    if smi:
      self.smilesRadio.setChecked(1)
      self._smilesLogged=0
      if self.constantUpdates:
        self.molUpdated()
        if chemdraw:
          self.toCDXButton.setEnabled(1)
  def smilesReturnPressed(self):
    self.molUpdated()
    self.logMols([(self.getMol(),self.getFormat())])
    self._smilesLogged=1
    
  def dragEnterEvent(self,evt):
    if QUriDrag.canDecode(evt):
      d = QStringList()
      QUriDrag.decodeLocalFiles(evt,d)
      ext = str(d.first()).split('.')[-1]
      if ext.lower() in self._supportedDropFormats:
        evt.accept(1)
        
  def dropEvent(self,evt):
    """ handle drop events

      **Notes**

        - if multiple files get dropped at once, we only look at the first

    """
    if QUriDrag.canDecode(evt):
      d = QStringList()
      if QUriDrag.decodeLocalFiles(evt,d):
        fileN = str(d.first())
      self.fileButtonClicked(fileN=fileN)

  def nextButtonClicked(self):
    done = 0
    while not done:
      try:
        self.nextMol()
      except IndexError:
        done=1
      else:
        for cb in self._nextMolCallbacks:
          cb(self.activeIdx())
        done = not self._processUpdatedMol(self.getMol(),self.getFormat())

  def previousButtonClicked(self):
    done = 0
    while not done:
      try:
        self.prevMol()
      except IndexError:
        done=1
      else:
        for cb in self._prevMolCallbacks:
          cb(self.activeIdx())
        done = not self._processUpdatedMol(self.getMol(),self.getFormat())

    

if __name__ == '__main__':
  from qtGui import Gui

  app,widg = Gui.Launcher(MolEditDisplayWidget,None,'MolEditDisplay',allowArrows=1,
                          logging=1)
  app.exec_loop()

