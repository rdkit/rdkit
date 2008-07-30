# $Id$
#
# Copyright (C) 2005,2006 Rational Discovery LLC
#  All Rights Reserved
#
from qt import *
from qtGui import qtUtils,GuiTable
from qtGui.GuiLib import MolCanvas
from Chem.Draw import MolDrawing
# FIX: this is a hack of the worst kind to get the boxes around
# labels in qt canvases smaller:
MolDrawing.MolDrawing.minLabelPadding=(0,-1)
import os,sys,types
import cPickle as pickle

import Chem
from Chem import ChemicalFeatures,ShowMols,rdDepictor,rdDistGeom
from qtGui.Search3D import SearchUtils,MolAlignmentList,LocalConfig
from qtGui.Search3D.PharmacophoreDefWidgets import FeatAssignmentWidget
import DistanceGeometry as DG

try:
  from utils import chemdraw
  reload(chemdraw)
except:
  import traceback
  traceback.print_exc()
  chemdraw=None

from qtGui.Search3D.forms.Search3DWidget import Search3DWidget as _Form
from qtGui.Search3D.PharmacophoreMixin import PharmacophoreMixin
from qtGui.Search3D.ExcludedVolumeMixin import ExcludedVolumeMixin

class Search3DWidget(_Form,PharmacophoreMixin,ExcludedVolumeMixin):
  def __init__(self,*args,**kwargs):
    _Form.__init__(self,*args,**kwargs)
    self._expertMode=True
    self._initCanvas()
    self._initContents()
    self._initViewer()
    self._initResultsPage()
    self._dir='.'
    
    self.slop = LocalConfig.distance3DSlop
    self.slop2D = LocalConfig.distance2DSlopFactor

    self.mol = None
    self.featAssignmentWidget=None
    self.distAssignmentWidget=None

    self.proteinName = ''
    self.proteinId=''

    self.excludedAtoms = {}
    self.useDirs=True
    self.addHs=False

    self.showProteinStateChanged()

  # --------------------------------------------------------------------
  #
  # Initialization/startup code:
  #
  # --------------------------------------------------------------------
  def _initCanvas(self):
    canvas = MolCanvas.MolCanvasView(self.loadPage)
    canvas.initCanvas((300,300))
    canvas.addRightHandler(lambda x,y:x.refresh())
    canvas.includeAtomNumbers=False
    self.refMolCanvas=canvas
    self.loadPage.layout().insertWidget(1,canvas)

    canvas = MolCanvas.MolCanvasView(self.pcophorePage)
    canvas.initCanvas((300,300))
    canvas.addRightHandler(lambda x,y:x.refresh())
    canvas.includeAtomNumbers=False
    self.pcophoreMolCanvas=canvas
    self.pcophorePage.layout().addWidget(canvas)

    
  def _initViewer(self):
    try:
      self.viewer3d = ShowMols.MolViewer(force=1)
    except:
      qtUtils.error('Could not contact PyMol application.\nPlease start PyMol and then rerun this application.')
      sys.exit(-1)
    if hasattr(self.viewer3d,'DisplayHBonds') and self.viewer3d.DisplayHBonds:
      self.supportHBonds=True
    else:
      self.supportHBonds=False
      self.showHBondsCheck.setEnabled(False)

  def _initResultsPage(self,sz=(300,300)):
    self.resultsSplitter = QSplitter(self.resultsPage)
    self.resultsSplitter.setOrientation(QSplitter.Vertical)
    #self.resultsPage.layout().insertWidget(2,self.resultsSplitter)
    self.resultsPage.layout().addWidget(self.resultsSplitter)

    self.molAlignList = MolAlignmentList.MolAlignmentList(self,self.resultsSplitter)


    canvas = MolCanvas.MolCanvasView(self.resultsSplitter)
    canvas.initCanvas(sz)
    
    canvas.addRightHandler(lambda x,y:x.refresh())
    canvas.includeAtomNumbers=False
    self.hitsMolCanvas=canvas

    if chemdraw:
      self.cdxGrabButton = QPushButton('Grab Molecule(s) from ChemDraw',self.resultsPage)
      #self.resultsPage.layout().insertWidget(3,self.cdxGrabButton)
      self.resultsPage.layout().addWidget(self.cdxGrabButton)
      QToolTip.add(self.cdxGrabButton,self.trUtf8("Select the active molecule(s) in ChemDraw and attempt to align them to the current pharmacophore."))
      self.connect(self.cdxGrabButton,SIGNAL('clicked()'),self.cdxGrab)
      self.cdxGrabButton.show()
    else:
      self.cdxGrabButton = QPushButton('Grab Smiles From Clipboard',self.resultsPage)
      self.resultsPage.layout().addWidget(self.cdxGrabButton)
      QToolTip.add(self.cdxGrabButton,self.trUtf8("Construct molecule(s) from SMILES string(s) on the clipboard and attempt to align them to the current pharmacophore."))
      self.connect(self.cdxGrabButton,SIGNAL('clicked()'),self.smilesGrab)
      self.cdxGrabButton.show()

  def _initContents(self):
    self.enablePages()
    if hasattr(self.topLevelWidget(),'enableDbWidget') and self.topLevelWidget().enableDbWidget:
      from qtGui.DbQueryWidgetImpl import insertQueryWidget

      layout10 = QHBoxLayout(None,0,2,"layout10")

      self.dbRadio = QRadioButton(self.inputTypeGroup,"dbRadio")
      self.dbRadio.setSizePolicy(QSizePolicy(0,0,0,0,self.dbRadio.sizePolicy().hasHeightForWidth()))
      layout10.addWidget(self.dbRadio)

      self.molDbWidgetBox = QGroupBox(self.inputTypeGroup,"molDbWidgetBox")
      layout10.addWidget(self.molDbWidgetBox)
      self.inputTypeGroup.layout().addLayout(layout10)

      self.molDbWidget = insertQueryWidget(self.molDbWidgetBox)
    else:
      self.molDbWidget = None
    self.loadMolsButton.setEnabled(True)
    self.addMolsButton.setEnabled(True)
    self.searchButton.setEnabled(False)
    self.enableHitButtons(False)

  def initFeatWidget(self):
    self.featAssignmentWidget = FeatAssignmentWidget(self.pcophorePage,self)
    self.pcophorePage.layout().insertWidget(2,self.featAssignmentWidget)
    self.featAssignmentWidget.show()
    #self.pcophorePage.layout().addStretch()

    self.grabDistancesButton=QPushButton(self.trUtf8('Grab Distances'),
                                         self.pcophorePage)
    self.connect(self.grabDistancesButton,
                 SIGNAL('clicked()'),
                 self.grabDistances)
    self.pcophorePage.layout().insertWidget(4,self.grabDistancesButton)
    self.grabDistancesButton.show()
    QToolTip.add(self.grabDistancesButton,
                 self.trUtf8("Initialize the pharmacophore with the distances between the current features."))
    #self.pcophorePage.layout().addStretch()

    self.initExcludeButtons()

  def updateCanvases(self):
    self.refMolCanvas.initCanvas()
    self.hitsMolCanvas.initCanvas()
    if hasattr(self,'pcophoreMolCanvas'):
      self.pcophoreMolCanvas.initCanvas()

  # --------------------------------------------------------------------
  #
  # Interface functionality:
  #
  # --------------------------------------------------------------------

  def enableHitButtons(self,state):
    self.refineSearchButton.setEnabled(state)
    self.clusterHitsButton.setEnabled(state)
    self.examineHitsButton.setEnabled(state)
    self.numPicksBox.setEnabled(state)

  def enablePages(self,loadStatus=True,pcophoreStatus=False,
                  detailsStatus=True,resultsStatus=True):
    if self.expertModeEnabled():
      self.tabWidget.setTabEnabled(self.pcophorePage,pcophoreStatus)
      self.tabWidget.setTabEnabled(self.loadPage,loadStatus)
      self.tabWidget.setTabEnabled(self.detailsPage,detailsStatus)
      self.tabWidget.setTabEnabled(self.resultsPage,resultsStatus)
    else:
      self.tabWidget.setTabEnabled(self.pcophorePage,True)
      self.tabWidget.setTabEnabled(self.resultsPage,True)

  def expertModeEnabled(self):
    return self._expertMode
    
  def enableExpertMode(self,state):
    expertWidgets=(self.fdefFileLabel,
                     self.fdefFileButton,
                     self.addAtomsButton,
                     self.fdefFilename,
                     self.centerViewCheck,
                     self.replaceCheck,
                     )
    
    self._expertMode=state
    if not state:
      for pg in (self.detailsPage,self.loadPage):
        if self.tabWidget.isTabEnabled(pg):
          pg._label = str(self.tabWidget.tabLabel(pg))
          pg._status = self.tabWidget.isTabEnabled(pg)
          self.tabWidget.removePage(pg)
      self.tabWidget.setTabEnabled(self.pcophorePage,True)
      self.tabWidget.setTabEnabled(self.resultsPage,True)

      for widget in expertWidgets:
        widget.hide()
    else:
      pg = self.loadPage
      if self.tabWidget.indexOf(pg)<0:
        self.tabWidget.insertTab(pg,pg._label,0)
      pg = self.detailsPage
      if self.tabWidget.indexOf(pg)<0:
        self.tabWidget.insertTab(pg,pg._label,2)

      for widget in expertWidgets:
        widget.show()
      if len(self.molAlignList.getMols()):
        self.tabWidget.setTabEnabled(self.resultsPage,True)
      if hasattr(self,'pcophore') and self.pcophore:
        self.tabWidget.setTabEnabled(self.pcophorePage,True)
        self.tabWidget.setTabEnabled(self.detailsPage,True)


    self.enableExcludeButtons(state)
    if hasattr(self,'activeFeats'):
      nFeats = len(self.activeFeats)
      self.featAssignmentWidget.enableRows(state,nFeats-1)
      self.distAssignmentWidget.enableRows(state,(nFeats*(nFeats-1))/2)
  def centerView(self):
    return self.centerViewCheck.isChecked()
  def replaceCurrent(self):
    return self.replaceCheck.isChecked()
  def showProtein(self):
    return self.showProteinCheck.isChecked()
  def showHBonds(self):
    return self.showHBondsCheck.isChecked()
  def showCollisions(self):
    return self.showCollisionsCheck.isChecked()


  # --------------------------------------------------------------------
  #
  # 3DViewer functionality:
  #
  # --------------------------------------------------------------------
  def clearAtomHighlights(self,where=LocalConfig.pymolTemplateName):
    if hasattr(self,'pcophoreMolCanvas'):
      self.pcophoreMolCanvas.clearHighlights()

  def highlightAtoms(self,indices,where=LocalConfig.pymolTemplateName,
                     extraHighlight=False):
    self.viewer3d.HighlightAtoms(indices,where,
                                 extraHighlight=extraHighlight)
    if hasattr(self,'pcophoreMolCanvas'):
      self.pcophoreMolCanvas.clearHighlights()
      self.pcophoreMolCanvas.highlightAtoms([x-1 for x in indices],
                                            highlightColor=(.7,.7,.7),
                                            highlightRadius=1.)

  def grabSelectedAtoms(self,tryPicksFirst=False,replacePicks=True):
    if not hasattr(self,'featureFactory') or not self.featureFactory:
      qtUtils.error('Please load an FDef file before grabbing or adding atoms')
      return
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if tryPicksFirst:
      selection = self.viewer3d.GetSelectedAtoms('pkset')
    else:
      selection = None
    if not selection:
      selection = self.viewer3d.GetSelectedAtoms()
    if self.expertModeEnabled() and (not selection or not len(selection)):
      QApplication.restoreOverrideCursor()
      qtUtils.error('no atoms selected in PyMol')
      return

    indices = []
    if selection:
      omolName = selection[0][0]
      for molName,idx in selection:
        if molName != omolName:
          QApplication.restoreOverrideCursor()
          qtUtils.error('atoms from multiple molecules selected')
          return
        else:
          indices.append(idx)

    if not self.expertModeEnabled():
      if len(indices)>LocalConfig.maxAddablePoints:
        qtUtils.error('You selected %d points, but only %d points\ncan be added to the pharmacophore.'%(len(indices),LocalConfig.maxAddablePoints))
        QApplication.restoreOverrideCursor()
        return

    if self.expertModeEnabled():
      if self.featAssignmentWidget is None:
        self.initFeatWidget()
      elif replacePicks:
        self.featAssignmentWidget.reinit()
      if replacePicks or not (hasattr(self,'activeIndices') and self.activeIndices):
        self.activeIndices = []

    else:
      if self.featAssignmentWidget is None:
        self.initFeatWidget()
      nLeft = self.featAssignmentWidget.removeEnabledRows()
      self.activeIndices = self.activeIndices[:nLeft]

    if self.updateFeatAssignmentWidget(indices):
      self.pcophoreMolCanvas.setMol(self.mol2d)

    QApplication.restoreOverrideCursor()

  def updateFeatAssignmentWidget(self,indices):
    if not hasattr(self,'featureFactory') or not self.featureFactory:
      return
    if self.featAssignmentWidget is None:
      self.initFeatWidget()

    feats = self.featureFactory.GetFeaturesForMol(self.mol)
    self.featMap = SearchUtils.GetFeatsPerAtom(feats)

    for entry in indices:
      if entry not in self.featAssignmentWidget.getAtomIds():
        try:
          nRows=self.featAssignmentWidget.addRow(entry,self.featMap[entry-1])
        except KeyError:
          self.highlightAtoms([entry],extraHighlight=True)
          qtUtils.error("No features found for atom %d.\nPlease check your selection."%entry)
          return False
        except:
          qtUtils.logger.debug("Couldn't add row",exc_info=True)
        else:
          self.activeIndices.append(entry)
      else:
        qtUtils.warning('A duplicate atom in the pharmacophore specification (atom %d) was ignored'%entry)

    return True
  def updateProtein(self):
    if not self.proteinName:
      self.proteinName = os.path.abspath(str(self.proteinFilename.text()))
    # FIX: handle removal of the protein:
    if self.proteinName:
      self.proteinId = self.viewer3d.LoadFile(self.proteinName,
                                              LocalConfig.pymolProteinName,
                                              showOnly=False)
      if not self.showProteinCheck.isChecked():
        self.viewer3d.HideObject(LocalConfig.pymolProteinName)
        self.viewer3d.Zoom(LocalConfig.pymolTemplateName)


  def displayHBonds(self,nm=LocalConfig.pymolHBondName,
                 mol=LocalConfig.pymolTemplateName,
                 protein=LocalConfig.pymolProteinName):
    self.viewer3d.DisplayHBonds(nm,mol,protein,
                                molSelText=LocalConfig.pymolMolSelText,
                                proteinSelText=LocalConfig.pymolProteinSelText)
                 
  def displayCollisions(self,nm=LocalConfig.pymolCollisionName,
                        mol=LocalConfig.pymolTemplateName,
                        protein=LocalConfig.pymolProteinName,
                        distCutoff=LocalConfig.pymolCollisionCutoff):
    self.viewer3d.DisplayCollisions(nm,mol,protein,distCutoff=distCutoff,
                                    molSelText=LocalConfig.pymolMolSelText,
                                    proteinSelText=LocalConfig.pymolProteinSelText)


  def cleanupViewer(self):
    res = self.molAlignList.getCheckedEmbeds(attr='getLabel')
    res.extend([LocalConfig.pymolProteinName,LocalConfig.pymolTemplateName,
                LocalConfig.refPcophoreName])
    self.viewer3d.DeleteAllExcept(res)

  # --------------------------------------------------------------------
  #
  # Chemical functionality:
  #
  # --------------------------------------------------------------------

  def loadMol(self,fileName=''):
    if not fileName:
      fileName = str(self.molFilename.text())
    try:
      inD = open(fileName,'r').read()
    except IOError:
      qtUtils.error('could not open file for reading:\n%s'%fileName)
      return

    try:
      self.mol = Chem.MolFromMolBlock(inD)
    except:
      qtUtils.error('could not construct molecule from file:\n%s'%fileName,
                    exc_info=True)
      return
    else:
      self.enablePages(True,True)
      # we no longer send the raw mol data to pymol, because
      # this was causing mismatched atom ids when the molecule
      # in question had hydrogens (RD Issue 359)
      #self.viewer3d.ShowMol(self.mol,
      #                      name=LocalConfig.pymolTemplateName,
      #                      molB=inD)
      self.viewer3d.ShowMol(self.mol,
                            name=LocalConfig.pymolTemplateName)

      self.viewer3d.SetDisplayStyle(LocalConfig.pymolTemplateName,
                                    'sticks')
      self.viewer3d.Zoom(LocalConfig.pymolTemplateName)

      self.mol2d = Chem.Mol(self.mol.ToBinary())
      Chem.Kekulize(self.mol2d)
      self.mol2d.RemoveAllConformers()
      rdDepictor.Compute2DCoords(self.mol2d)
      #print 'load set'
      self.refMolCanvas.setMol(self.mol2d)

  def loadFDefs(self,fileName='',noWarning=False,noErrorMessages=False):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if not fileName:
      fileName = str(self.fdefFilename.text())
    try:
      inD = open(fileName,'r').read()
    except IOError:
      QApplication.restoreOverrideCursor()
      if not noErrorMessages:
        qtUtils.error('could not open file for reading:\n%s'%fileName)
      return
    try:
      self.featureFactory = ChemicalFeatures.BuildFeatureFactoryFromString(inD)
    except:
      QApplication.restoreOverrideCursor()
      if not noErrorMessages:
        qtUtils.error('could not load feature definitions from file:\n%s'%fileName,
                      exc_info=True)
      self.featureFactory = None
      return
    QApplication.restoreOverrideCursor()

    if not noWarning and hasattr(self,'pcophore') and self.pcophore:
      confirm = QMessageBox.warning(self,
                                    self.trUtf8("Discard Results?"),
                                    self.trUtf8("Loading a new FDef file will discard your current results. Continue?"),
                                    1,2)
      if confirm != 1:
        return

    if self.featAssignmentWidget is not None:
      self.featAssignmentWidget.reinit()
    self.pcophore=None
    if hasattr(self,'mol2d') and self.mol2d:
      self.pcophoreMolCanvas.setMol(self.mol2d)
    if self.distAssignmentWidget:
      self.distAssignmentWidget.reinit()
    if hasattr(self,'molAlignList') and self.molAlignList:
      self.molAlignList.clear()

  def addTargets(self):
    if self.sdFileRadio.isChecked():
      self.addTargetsFromFile()
    elif hasattr(self,'dbRadio') and self.dbRadio.isChecked():
      self.addTargetsFromDb()

  def loadTargets(self):
    if hasattr(self,'targets'):
      saveTgts = self.targets
    else:
      saveTgts = []
    self.targets = []
    if self.sdFileRadio.isChecked():
      self.addTargetsFromFile()
    elif hasattr(self,'dbRadio') and self.dbRadio.isChecked():
      self.addTargetsFromDb()
    # recover from error:
    if not self.targets:
      self.targets = saveTgts

  def addTargetsFromDb(self,conn=None):
    from Chem.Suppliers.DbMolSupplier import ForwardDbMolSupplier as DbMolSupplier
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    rawD = self.molDbWidget.getData(randomAccess=0)
    try:
      supplier = DbMolSupplier(rawD,
                               nameCol='ID')
    except ValueError:
      qtUtils.warning('Problems encountered loading compounds from database.\nPlease check that you have selected a table containing molecules.',
                      exc_info=True)
      QApplication.restoreOverrideCursor()
      return
      
    count = self.molDbWidget.getDataCount()
    self.addTargetsFromSupplier(supplier,count=count)
    QApplication.restoreOverrideCursor()
  
  def addTargetsFromFile(self,fileName=''):
    if not fileName:
      fileName = str(self.sdFilename.text())
    try:
      inD = open(fileName,'r').read()
    except IOError:
      qtUtils.error('could not open file for reading:\n%s'%fileName)
      return
    try:
      supplier = Chem.SupplierFromFilename(fileName,sanitize=False)
    except ValueError:
      qtUtils.error('Unrecognized file type for file:\n%s'%fileName)
      return
    self.addTargetsFromSupplier(supplier)

  def addTargetsFromSupplier(self,supplier,count=-1):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if count < 0:
      count = len(supplier)
    dlg = qtUtils.ProgressDialogHolder('Loading Molecules',count)
    if hasattr(self,'targets') and self.targets:
      targets = self.targets
    else:
      targets=[]
    try:
      SearchUtils.InitializeTargets(supplier,progressCallback=dlg,
                                    res=targets,addHs=self.addHs)
    except qtUtils.ProgressStop:
      qtUtils.warning("Aborted")
    except:
      qtUtils.logger.debug("Couldn't initialize targets",exc_info=True)
    if targets:
      self.targets=targets
      nMols = len(self.targets)
      self.numMolsLabel.setText('Number of Molecules: %d'%nMols)
      self.searchButton.setEnabled(True)

    QApplication.restoreOverrideCursor()

  def searchMols(self,pcophore=None,excludes=[],targets=None):
    if not targets and hasattr(self,'targets'):
      targets = self.targets
    if not targets:
      return
    if not pcophore:
      pcophore = self.buildPharmacophore()
    if not pcophore:
      return
    if not excludes:
      excludes = self.buildExcludedVolumes()
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    count = len(targets)
    dlg = qtUtils.ProgressDialogHolder('Searching Molecules',count)
    try:
      hits = []
      SearchUtils.SearchTargets(targets,self.featureFactory,
                                pcophore,excludedVolumes=excludes,
                                useDirs=self.useDirs,
                                progressCallback=dlg,
                                res=hits)
    except qtUtils.ProgressStop:
      QApplication.restoreOverrideCursor()
      qtUtils.warning("Aborted")
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.error('error encountered while searching',exc_info=True)
      return
    if hits:
      self.hits=hits
      nHits = len(self.hits)
      self.numHitsLabel.setText('Number of Hits: %d'%nHits)
      if nHits:
        self.enableHitButtons(True)
        self.numPicksBox.setMinValue(1)
        self.numPicksBox.setMaxValue(nHits)
        currVal = min(10,nHits/2)
        self.numPicksBox.setValue(currVal)
      else:
        self.enableHitButtons(False)
    else:
      qtUtils.warning('No matching molecules found.')
    self.enablePages(True,True,True)
    QApplication.restoreOverrideCursor()

  def transferHits(self,clearExisting=True):
    if not hasattr(self,'hits') or not self.hits:
      qtUtils.error("No hits to transfer")
      return
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    nHits = len(self.hits)
    #self.numHitsLabel2.setText('Number of Hits: %d'%nHits)

    if clearExisting:
      self.molAlignList.clear()

    dlg = qtUtils.ProgressDialogHolder('Transferring Results',nHits)
    try:
      for i in range(nHits):
        mol = self.hits[i][0]
        bounds = self.hits[i][-1]
        idx = self.molAlignList.childCount()+1
        if mol.HasProp('_Name'):
          label = mol.GetProp('_Name')
        else:
          label = ''
        if not label:
          label = '%s-%d'%(LocalConfig.resultsMoleculeNameBase,idx)
          mol.SetProp('_Name',label)
        item = MolAlignmentList.MolItem(mol,bounds,label,idx,self.molAlignList)
        item.setVisible(1)
        if not i:
          self.molAlignList.ensureItemVisible(item)
          self.molAlignList.clearSelection()
          self.molAlignList.setSelected(item,True)
        dlg(i)
    except qtUtils.ProgressStop:
      QApplication.restoreOverrideCursor()
      qtUtils.warning("Aborted")
    else:
      self.molAlignList.show()
    self.enablePages(True,True,True,True)
    QApplication.restoreOverrideCursor()

  def getSearchResults(self,checkedOnly=False):
    if checkedOnly:
      res = self.molAlignList.getCheckedEmbeds()
    else:
      res = self.molAlignList.getMols()
    return res

  def convertHits(self):
    if not hasattr(self,'hits') or not self.hits:
      return
    self.targets = []
    for tpl in self.hits:
      hit = tpl[0]
      bm = rdDistGeom.GetMoleculeBoundsMatrix(hit)
      if DG.DoTriangleSmoothing(bm):
        self.targets.append((hit,bm))
      else:
        qtUtils.logger.info('problems encountered smoothing bounds matrix for molecule %s'%hit.GetProp('_Name'))
        
    self.hits = []
    self.numMolsLabel.setText('Number of Molecules: %d'%len(self.targets))
    self.enableHitButtons(False)
    self.searchButton.setEnabled(True)
    
  def clusterHits(self):
    if not hasattr(self,'hits') or not self.hits:
      return
    nToGrab = self.numPicksBox.value()
    if nToGrab >= len(self.hits):
      return

    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    nHits= len(self.hits)
    nElems = nHits*(nHits-1)/2
    dlg = qtUtils.ProgressDialogHolder('Generating Distance Matrix',nElems)
    try:
      picks = SearchUtils.DiversityPick(self.hits,nToGrab,progressCallback=dlg)
    except qtUtils.ProgressStop:
      QApplication.restoreOverrideCursor()
      qtUtils.warning("Aborted")
    except:
      qtUtils.logger.debug("Couldn't do diversity pick",exc_info=True)
    else:
      self.hits = picks
      nHits = len(self.hits)
      self.numHitsLabel.setText('Number of Hits: %d'%nHits)
      self.numPicksBox.setMaxValue(nHits)
      currVal = min(10,nHits/2)
      self.numPicksBox.setValue(currVal)
    QApplication.restoreOverrideCursor()


  # --------------------------------------------------------------------
  #
  # Signals and Slots
  #
  # --------------------------------------------------------------------
  def molFilenameChanged(self):
    fName = str(self.molFilename.text())
    self.loadMol(fileName=fName)

  def molFileButtonClicked(self):
    fileN = str(QFileDialog.getOpenFileName(self._dir,
                                            'MDL Mol files (*.mol)'))
    if fileN:
      self.molFilename.setText(fileN)
      self._dir = os.path.dirname(fileN)
      self.loadMol(fileName=fileN)


  def fdefFilenameChanged(self):
    fName = str(self.fdefFilename.text())
    self.loadFDefs(fileName=fName)

  def fdefFileButtonClicked(self):
    fileN = str(QFileDialog.getOpenFileName(self._dir,
                                            'Feature definition files (*.fdef)'))
    if fileN:
      self.fdefFilename.setText(fileN)
      self._dir = os.path.dirname(fileN)
      self.loadFDefs(fileName=fileN)


  def grabAtomsButtonClicked(self):
    self.grabSelectedAtoms()
  def addAtomsButtonClicked(self):
    self.grabSelectedAtoms(replacePicks=False)



  def sdFilenameChanged(self):
    fName = str(self.sdFilename.text())
    #self.loadMolsButton.setEnabled(True)
    #self.loadTargetsFromFile(fileName=fName)

  def sdFileButtonClicked(self):
    fileN = str(QFileDialog.getOpenFileName(self._dir,
                                            'SD files (*.sdf)'))
    if fileN:
      self._dir = os.path.dirname(fileN)
      self.loadMolsButton.setEnabled(True)
      self.addMolsButton.setEnabled(True)
      self.sdFilename.setText(fileN)
      #self.addTargetsFromFile(fileName=fileN)


  def loadSdFileAction(self,id):
    fileN = str(QFileDialog.getOpenFileName(self._dir,
                                            'SD files (*.sdf)'))
    if fileN:
      self.targets=[]
      self._dir = os.path.dirname(fileN)
      self.loadMolsButton.setEnabled(True)
      self.addMolsButton.setEnabled(True)
      self.sdFilename.setText(fileN)
      self.addTargetsFromFile(fileName=fileN)
      self.searchMols()
      self.transferHits(clearExisting=False)
  def loadMolsButtonClicked(self):
    self.loadTargets()
  def addMolsButtonClicked(self):
    self.addTargets()

  def searchButtonClicked(self):
    self.searchMols()

  def refineSearchButtonClicked(self):
    if not hasattr(self,'hits') or not self.hits:
      return
    self.convertHits()

  def clusterHitsButtonClicked(self):
    if not hasattr(self,'hits') or not self.hits:
      return
    self.clusterHits()

  def examineHitsClicked(self):
    self.transferHits()
    self.tabWidget.setCurrentPage(3)

  def addExclusionClicked(self):
    self.addExcludedAtoms()

  def resetExclusionClicked(self):
    self.resetExcludedAtoms()

  def proteinFilenameChanged(self):
    fName = str(self.proteinFilename.text())
    self.proteinName = fName
    self.proteinId = ''
    self.updateProtein()
  def proteinFileButtonClicked(self):
    fileN = str(QFileDialog.getOpenFileName(self._dir,
                                            "PDB files (*.pdb)"))
    if fileN:
      self._dir = os.path.dirname(fileN)
      self.proteinFilename.setText(fileN)
      self.proteinName = fileN
      self.proteinId = ''
      self.updateProtein()

  def neighborStateChanged(self):
    if self.neighboringSurfaceCheck.isChecked():
      self.viewer3d.SetDisplayStyle(LocalConfig.pymolProteinName,'')
      self.viewer3d.SelectProteinNeighborhood(LocalConfig.pymolTemplateName,
                                              LocalConfig.pymolProteinName,
                                              distance=LocalConfig.neighborhoodRadius,
                                              name='neighborhood',showSurface=True)
    else:
      self.viewer3d.SetDisplayStyle('neighborhood','')
      self.viewer3d.HideObject('neighborhood')
      if self.showProtein():
        self.viewer3d.SetDisplayStyle(LocalConfig.pymolProteinName,
                                      'lines')

  def showProteinStateChanged(self):
    if self.showProtein():
      self.viewer3d.DisplayObject(LocalConfig.pymolProteinName)
      self.viewer3d.SetDisplayStyle(LocalConfig.pymolProteinName,
                                    'lines')
      self.neighboringSurfaceCheck.setEnabled(True)
      if self.supportHBonds and (not LocalConfig.hBondsAlwaysEnabled):
        self.showHBondsCheck.setEnabled(True)
      if self.neighboringSurfaceCheck.isChecked():
        self.neighborStateChanged()
    else:
      self.neighboringSurfaceCheck.setEnabled(False)
      if self.supportHBonds and (not LocalConfig.hBondsAlwaysEnabled):
        self.showHBondsCheck.setEnabled(False)
      self.viewer3d.HideObject(LocalConfig.pymolProteinName)


  def showHBondsStateChanged(self):
    if self.showHBonds():
      itm = self.molAlignList.currentItem()
      if itm and isinstance(itm,MolAlignmentList.EmbedItem) and \
         itm.dispLabel:
        self.displayHBonds(mol=itm.dispLabel)
    else:
      self.viewer3d.HideObject(LocalConfig.pymolHBondName)
      #self.viewer3d.server.do('disable %s'%LocalConfig.pymolHBondName)

  def showCollisionsStateChanged(self):
    if self.showCollisions():
      itm = self.molAlignList.currentItem()
      if itm and isinstance(itm,MolAlignmentList.EmbedItem) and \
         itm.dispLabel:
        self.displayCollisions(mol=itm.dispLabel)
    else:
      self.viewer3d.HideObject(LocalConfig.pymolCollisionName)
      #self.viewer3d.server.do('disable %s'%LocalConfig.pymolHBondName)



  def cdxGrab(self,id=-1,doSearch=LocalConfig.searchOnCDXGrab):
    if not hasattr(self,'activeFeats') or not self.activeFeats:
      qtUtils.error('Please define pharmacophore features and distances before grabbing molecules')
      return
    if not hasattr(self,'distAssignmentWidget') or not self.distAssignmentWidget or \
       not hasattr(self.distAssignmentWidget,'rows') or \
       not self.distAssignmentWidget.rows:
      qtUtils.error('Please define pharmacophore distances before grabbing molecules')
      return
      
    if not chemdraw:
      return
    elif not hasattr(chemdraw,'cdApp'):
      qtUtils.warning('Problems were encountered connecting to Chemdraw.')
      return
    try:
      if chemdraw.cdApp is None:
        chemdraw.StartChemDraw()
      d = str(chemdraw.CDXGrab(outFormat='chemical/x-daylight-smiles'))
    except:
      qtUtils.warning('Could not grab molecule from Chemdraw',exc_info=True)
      return
    if not d:
      qtUtils.logger.debug('no data available from Chemdraw')
      return

    self.smilesGrab(id=id,doSearch=doSearch,data=d)
      
  def smilesGrab(self,id=-1,doSearch=LocalConfig.searchOnCDXGrab,data=None):
    if not hasattr(self,'activeFeats') or not self.activeFeats:
      qtUtils.error('Please define pharmacophore features and distances before grabbing molecules')
      return
    if not hasattr(self,'distAssignmentWidget') or not self.distAssignmentWidget or \
       not hasattr(self.distAssignmentWidget,'rows') or \
       not self.distAssignmentWidget.rows:
      qtUtils.error('Please define pharmacophore distances before grabbing molecules')
      return
      
    if data is None:
      data = str(QApplication.clipboard().text())
      if not data:
        qtUtils.warning('No text on clipboard')
        return
    smis = data.split('.')
    mols = []
    numProblems=0
    for smi in smis:
      try:
        m = Chem.MolFromSmiles(smi)
        mols.append(m)
      except:
        qtUtils.logger.debug("Couldn't build molecule from smiles: %s"%smi,
                             exc_info=True)
        numProblems +=1 
    if numProblems >0:
      quantity = '%d molecule'%numProblems
      if numProblems > 1:
        quantity += 's'
      msg = 'Problems were encountered building %s.\nPlease check the data to ensure that they are correct.\n'%quantity
      if len(mols):
        msg += 'The other molecules were successfully added.'
      qtUtils.warning(msg)
    if not len(mols):
      return
    targets =[]
    SearchUtils.InitializeTargets(mols,res=targets,addHs=self.addHs)
    pcophore = self.buildPharmacophore()
    excludes = self.buildExcludedVolumes()
    if doSearch:
      hits = []
      SearchUtils.SearchTargets(targets,self.featureFactory,
                                pcophore,excludedVolumes=excludes,
                                useDirs=self.useDirs,
                                res=hits)
      if not hits:
        qtUtils.information("no hits found among the new molecules",
                            caption="Search Results")
    else:
      hits = targets
    if not hits:
      return
    else:
      qtUtils.information("%d new molecules found"%len(hits),
                          caption="Search Results")
      for i,hit in enumerate(hits):
        mol = hit[0]
        bounds = hit[-1]
        if mol.HasProp('_Name'):
          label = mol.GetProp('_Name')
        else:
          label = ''
        idx = self.molAlignList.childCount()+1
        if not label:
          label = '%s-%d'%(LocalConfig.resultsMoleculeNameBase,idx)
          mol.SetProp('_Name',label)
        item = MolAlignmentList.MolItem(mol,bounds,label,idx,self.molAlignList)
        item.setVisible(1)
        if not i:
          self.molAlignList.ensureItemVisible(item)
          self.molAlignList.clearSelection()
          self.molAlignList.setSelected(item,True)
      
  # --------------------------------------------------------------------
  #
  # Serialization (pickling) support:
  #
  # --------------------------------------------------------------------
  def toString(self):
    from qtGui.Search3D import Persist
    res = {}

    # "globals":
    res['configModule'] = Persist.ModuleVarsToDict(LocalConfig)

    # things from the GUI:
    res['molFilename'] = str(self.molFilename.text())
    res['molFileData'] = open(str(self.molFilename.text()),'r').read()
    res['proteinFilename'] = str(self.proteinFilename.text())
    if res['proteinFilename']:
      res['proteinFileData'] = open(str(self.proteinFilename.text()),'r').read()
    else:
      res['proteinFileData'] = ''
    res['fdefFilename'] = str(self.fdefFilename.text())
    res['fdefData'] = open(str(self.fdefFilename.text()),'r').read()
    
    # calculated objects:
    if hasattr(self,'activeIndices'):
      res['activeIndices'] = pickle.dumps(self.activeIndices,True)
      res['activeFeats'] = pickle.dumps([x.GetFamily() for x in self.activeFeats],True)
      res['numEnabledRows'] = self.featAssignmentWidget.numEnabledRows()
    else:
      res['activeIndices'] = pickle.dumps([],True)
      res['activeFeats'] = pickle.dumps([],True)
      res['numEnabledRows'] = 0
    
    if not hasattr(self,'pcophore'):
      self.buildPharmacophore()
    res['pcophore'] = pickle.dumps(self.pcophore,True)
    if hasattr(self,'excludedVols'):
      res['excludedVols'] = pickle.dumps(self.excludedVols,True)
    else:
      res['excludedVols'] = ''
      
    

    # molecules and bounds matrices:
    items = self.molAlignList.getAllInstancesOf(MolAlignmentList.MolItem)
    mols = []
    bounds = []
    children = []
    for item in items:
      mol = item.mol
      name = mol.GetProp('_Name')
      mols.append((name,mol))
      bounds.append(item.bounds)
      children.append((item._energy,item.getChildState()))
    res['molList'] = mols
    res['boundsMatrixList']= bounds
    res['molListChildren'] = children

    return pickle.dumps(res,True)

  def fromString(self,pkl):
    import time
    modeSave = self._expertMode
    self._expertMode=True
    global tempHandler
    from qtGui.Search3D import Persist
    
    tempHandler = Persist.TempFileHandler()
    obj = pickle.loads(pkl)
    if type(obj) is not types.DictType:
      raise ValueError,'bad pickle'

    Persist.DictToModuleVars(obj['configModule'],LocalConfig)
    molFilename = tempHandler.get('.mol')
    open(molFilename,'w+').write(obj['molFileData'])
    self.molFilename.setText(molFilename)
    self.loadMol()
    
    proteinD = obj.get('proteinFileData','')
    if proteinD:
      proteinFilename = tempHandler.get('.pdb')
      outF = file(proteinFilename,'w+')
      outF.write(proteinD)
      outF.close()
      outF = None
      time.sleep(0.5)
      self.proteinFilename.setText(proteinFilename)
      self.proteinName=proteinFilename
      self.updateProtein()

    fdefFilename = tempHandler.get('.fdef')
    open(fdefFilename,'w+').write(obj['fdefData'])
    self.fdefFilename.setText(fdefFilename)
    self.loadFDefs(noWarning=True,noErrorMessages=True)
    if not hasattr(self,'featureFactory') or self.featureFactory is None:
      err = 'Without feature definitions, this session is unuseable.\n'
      if hasattr(sys,'frozen') and sys.frozen==1:
        err += 'Please send the log file to your support staff.'
      else:
        err += 'Exiting.'
      qtUtils.error(err)
      sys.exit(-1)
      
      

    indices=pickle.loads(obj['activeIndices'])
    self.highlightAtoms(indices)
    self.activeIndices=[]
    if self.featAssignmentWidget is not None:
      self.featAssignmentWidget.reinit()

    self.updateFeatAssignmentWidget(indices)
    idx = self.tabWidget.indexOf(self.pcophorePage)
    self.tabWidget.setCurrentPage(idx)
    self.pcophoreMolCanvas.setMol(self.mol2d)

    activeFeats = pickle.loads(obj['activeFeats'])
    self.pcophore = pickle.loads(obj['pcophore'])
    if self.pcophore:
      if len(activeFeats) != len(self.pcophore.getFeatures()):
        raise ValueError,'feature length mismatch'
      for i,featName in enumerate(activeFeats):
        box = self.featAssignmentWidget.rows[i][0]
        for j in range(box.count()):
          box.setCurrentItem(j)
          if str(box.currentText())==featName:
            break
      
    volPkl = obj.get('excludedVols','')
    if volPkl:
      self.excludedVols = pickle.loads(volPkl)
    else:
      self.excludedVols = None
    if activeFeats:
      self.grabDistances()
      offset=1
      for i in range(len(self.activeFeats)):
        for j in range(i+1,len(self.activeFeats)):
          row = self.distAssignmentWidget.rows[offset]
          offset += 1
          row[0].setText(str(self.pcophore.getLowerBound(i,j)))
          row[1].setText(str(self.pcophore.getUpperBound(i,j)))
          row[2].setValue(int(self.pcophore.getUpperBound2D(i,j)))


    mols = obj.get('molList',[])
    bounds = obj.get('boundsMatrixList',[])
    children = obj.get('molListChildren',[])
    for idx in range(len(mols)):
      nm,mol = mols[idx]
      mol.SetProp('_Name',nm)
      bm = bounds[idx]
      item = MolAlignmentList.MolItem(mol,bm,nm,idx,self.molAlignList)
      if children and children[idx]:
        energy,state= children[idx]
        item._energy=energy
        item.restoreChildState(state)
      if not idx:
        self.molAlignList.ensureItemVisible(item)
        self.molAlignList.clearSelection()
        self.molAlignList.setSelected(item,True)
    self.enableExpertMode(modeSave)

    firstEnabled = len(activeFeats)-obj.get('numEnabledRows',len(activeFeats))
    if firstEnabled>0:
      for i in range(firstEnabled,len(activeFeats)):
        self.featAssignmentWidget.enableRow(i,True)

    self.viewer3d.Zoom(LocalConfig.pymolTemplateName)
    self.displayPharmacophore()

  def __getstate__(self):
    return self.toString()

  # --------------------------------------------------------------------
  #
  # Debug code
  #
  # --------------------------------------------------------------------
  def cdk2Init(self,*args):
    if self.molDbWidget:
      self.molDbWidget.nameBox.setText("::ChemCatalogs")
      self.molDbWidget.tableCombo.insertItem("SmallCompoundsBin",0)
      self.molDbWidget.tableCombo.setCurrentItem(0)
      self.molDbWidget.enableQueries()

    self.molFilename.setText('testData/1oir-xtal.mol')
    self.proteinFilename.setText('testData/1OIR-nowater.pdb')
    self.loadMol()
    self.updateProtein()
    self.tabWidget.setCurrentPage(1)
    #self.viewer3d.server.do('zoom %s'%LocalConfig.pymolTemplateName)
    self.fdefFilename.setText('testData/SimpleFeatures.fdef')
    #self.fdefFilename.setText('testData/BaseFeatures.fdef')
    self.loadFDefs()

    self.sdFilename.setText('testData/cdk2-syn-clip100.sdf')
    self.loadMolsButton.setEnabled(True)
    self.addMolsButton.setEnabled(True)


    self.tabWidget.setCurrentPage(2)
    # CDK2:
    self.highlightAtoms((14,17,22))
    #self.highlightAtoms((14,17,22,29))
    self.grabSelectedAtoms()

    self.featAssignmentWidget.rows[0][0].setCurrentItem(1)
    self.featAssignmentWidget.rows[1][0].setCurrentItem(1)
    self.featAssignmentWidget.rows[2][0].setCurrentItem(1)
    if len(self.featAssignmentWidget.rows)>3:
      self.featAssignmentWidget.rows[3][0].setCurrentItem(2)
    self.grabDistances()
    self.distAssignmentWidget.rows[2][2].setValue(6)
    self.tabWidget.setCurrentPage(1)


  def sarsInit(self,*args):
    if self.molDbWidget:
      self.molDbWidget.nameBox.setText("::ChemCatalogs")
      self.molDbWidget.tableCombo.insertItem("SmallCompoundsBin",0)
      self.molDbWidget.tableCombo.setCurrentItem(0)
      self.molDbWidget.enableQueries()

    self.molFilename.setText('testData/1P9U_ligand.mol')
    self.proteinFilename.setText('./testData/1P9U-ChainB.pdb')
    self.loadMol()
    self.updateProtein()
    self.viewer3d.server.do('zoom %s'%LocalUtils.pymolTemplateName)
    self.fdefFilename.setText('testData/SimpleFeaturesfdef')
    self.fdefFilename.setText('testData/BaseFeatures.fdef')
    self.loadFDefs()

    self.sdFilename.setText('testData/cdk2-syn-clip100.sdf')
    self.loadMolsButton.setEnabled(True)
    self.addMolsButton.setEnabled(True)

    self.tabWidget.setCurrentPage(2)

    # SARS:
    self.highlightAtoms((14,15,16,22))
    
    self.grabSelectedAtoms()
    self.featAssignmentWidget.rows[0][0].setCurrentItem(1)
    self.featAssignmentWidget.rows[1][0].setCurrentItem(1)
    self.featAssignmentWidget.rows[2][0].setCurrentItem(1)
    if len(self.featAssignmentWidget.rows)>3:
      self.featAssignmentWidget.rows[3][0].setCurrentItem(1)
    self.grabDistances()
    self.distAssignmentWidget.rows[2][2].setValue(6)
    self.tabWidget.setCurrentPage(1)
    self.viewer3d.server.do('enable Protein')
    
  def debugInit(self):
    win = self.topLevelWidget()
    win._actionMenu.insertItem(self.trUtf8('&CDK2 Init'),self.cdk2Init)
    win._actionMenu.insertItem(self.trUtf8('&SARS Init'),self.sarsInit)
    if 1:
      txt = open('save1.p3d','rb').read()
      self.fromString(txt)
    else:
      self.cdk2Init()
      open('cdk2_state.pkl','wb+').write(self.toString())
    self.tabWidget.setCurrentPage(1)
    return
    self.tabWidget.setCurrentPage(2)
    self.loadTargets()
    self.searchMols()

        
    #return
    self.transferHits()
    self.tabWidget.setCurrentPage(3)

    

