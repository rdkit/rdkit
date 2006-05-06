# $Id: MolAlignmentList.py 5084 2006-03-15 15:01:25Z glandrum $
#
# Copyright (C) 2005,2006 Rational Discovery LLC
#  All Rights Reserved
#
from Numeric import *
from qt import *
from qtGui import qtUtils
import cPickle
import Chem
from Chem import ChemicalForceFields,rdDepictor,rdDistGeom
from Chem.Pharm3D import EmbedLib,Pharmacophore,ExcludedVolume
EmbedLib.logger = qtUtils.logger
from qtGui.Search3D import SearchUtils,LocalConfig
from qtGui.Search3D.AtomCorrespondenceDialog import AtomCorrespondenceDialog

try:
  from utils import chemdraw
  reload(chemdraw)
except ImportError:
  hasCDX=0
else:
  try:
    tmp = chemdraw.CDXDisplay
  except AttributeError:
    hasCDX=0
  else:
    hasCDX=1
    import win32gui,time
cdApp=None

# Qt doesn't normally provide custom tooltips on QListViewItems,
# so we need to create a special class for that:
class MolAlignmentListToolTip(QToolTip):
  def __init__(self,parent):
    QToolTip.__init__(self,parent.viewport())
    self.listView = parent
    self._tipDict = {}
  def add(self,widg,mesg):
    self._tipDict[widg]=mesg
  def remove(self,widg):
    if self._tipDict.has_key(widg):
      del self._tipDict[widg]
      
  def maybeTip(self,pt):
    if not self.listView:
      return
    itm = self.listView.itemAt(pt)
    if not itm:
      return

    itmText = self._tipDict.get(itm,'')
    if not itmText:
      return

    itmRect = self.listView.itemRect(itm)
    if not itmRect.isValid():
      return

    col = self.listView.header().sectionAt(pt.x())
    # general purpose test:
    if col==-1:
      return

    # particular for us:
    if col != 0:
      return
    

    headerRect=self.listView.header().sectionRect(col)
    if not headerRect.isValid():
      return

    cellRect = QRect(headerRect.left(),itmRect.top(),
                     headerRect.width(),itmRect.height())
    self.tip(cellRect,itmText)
    
    


class MolAlignmentListItem(QListViewItem):
  """ general-purpose functionality for entries in the MolAlignmentList

  We implement compare() and text()

  """
  def __init__(self,label,keyVal,*args):
    QListViewItem.__init__(self,*args)
    self.label = label.replace(' ','_').replace('\t','_')
    self.keyVal = keyVal
    self.setRenameEnabled(0,True)
    self.setText(0,label)
  def compare(self,other,col,ascending):
    if not isinstance(other,MolAlignmentListItem) or col!=0:
      return 0
    else:
      return cmp(self.keyVal,other.keyVal)

  def text(self,col):
    if col==0:
      return QListViewItem.text(self,0)
    else:
      return '...'

  def okRename(self,col):
    QListViewItem.okRename(self,col)
    if col==0:
      self.label = str(self.text(0))

  def getName(self):
    """ looks up the name of the molecule (which may have changed due to
    tree edits) and returns it
    """
    p = self.parent()
    op = self
    while p:
      op = p
      p = self.parent()
    nm = str(op.text(0))
    return nm
  
  def getMol(self):
    """ looks up the name of the molecule (which may have changed due to
    tree edits), sets the name property, and returns the molecule
    """
    self.mol.SetProp('_Name',self.getName())
    return self.mol
  def getMol2d(self):
    """ looks up the name of the molecule (which may have changed due to
    tree edits), sets the name property, and returns the molecule
    """
    self.mol2d.SetProp('_Name',self.getName())
    return self.mol2d


class EmbedItem(QCheckListItem):
  """ A class for storing actual embeddings 

  Here our compare() and text() methods are more sophisticated

  These are leaf nodes in the tree, they cannot be further expanded.

  """
  def __init__(self,mapping,mol,label,keyVal,exVols,parent,**kwargs):
    """ 

    """
    QCheckListItem.__init__(self,parent,label,QCheckListItem.CheckBox)
    self.label = label
    self.keyVal = keyVal
    self.mol = mol
    self.mapping = mapping
    self.exVols = exVols
    self.aligned=False
    self.setExpandable(False)
    self.index = kwargs.get('index',int(keyVal))
    if not kwargs.get('preAligned',False):
      self.initChemistry()
    self.dispLabel=''
    self.contextMenu=None

  def initChemistry(self):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    owner = self.listView().owner
    ref = owner.mol
    refFeats = owner.activeFeats
    probeFeats = SearchUtils.GetFeaturesForAtomMatch(self.mol,self.mapping,
                                                     owner.pcophore,owner.featureFactory)
    if len(probeFeats)!= len(refFeats):
      qtUtils.error('feature length mismatch')
      return

    rms,volPs,shapeScore,tform = SearchUtils.AlignMatchToReference(self.mol,probeFeats,
                                                                   ref,refFeats,
                                                                   exVols=self.exVols,
                                                                   useDirs=owner.useDirs)
    self.rms = rms
    self.shapeScore=shapeScore
    self.mol.SetProp('AlignmentRMS',str(rms))
    self.mol.SetProp('AlignmentShapeScore',str(shapeScore))
    
    self.volPs= volPs
    self.tform = tform
    self.aligned=True
    QApplication.restoreOverrideCursor()

  def contextMenuRequested(self,itm,pos,col):
    if self!=itm:
      return
    if not self.contextMenu:
      self.contextMenu = QPopupMenu(self.listView())
      self.contextMenu.insertItem('&Refine Alignment',self.refineAlignment)
    menu = self.contextMenu
    if menu:
      res = menu.exec_loop(pos)

  def refineAlignment(self):
    if not hasattr(self,'atomCorrespDialog') or not self.atomCorrespDialog:
      self.atomCorrespDialog = AtomCorrespondenceDialog(self,
                                                        self.listView().owner.mol,
                                                        #parent=self.listView().topLevelWidget())
                                                        parent=None)
    else:
      self.atomCorrespDialog.setProbeMol(self.getMol())
      self.atomCorrespDialog.initCorrespondence()
    self.atomCorrespDialog.show()
      
    
  def activate(self):
    if not self.isSelected():
      return
    QListViewItem.activate(self)
    if not self.aligned:
      self.initChemistry()
      
    owner = self.listView().owner
    # update the 2D drawing:
    self.parent().drawMol()
    qtUtils.logger.debug('RMS: %f'%self.rms)
    label = '%s-%d-%d'%(self.parent().parent().label,self.parent().keyVal,self.index)
    self.dispLabel=label
    owner.viewer3d.SetDisplayUpdate(False)
    if owner.replaceCurrent():
      owner.viewer3d.HideAll()
      owner.viewer3d.ShowMol(self.mol,name=label,showOnly=False,
                             zoom=owner.centerView())
      show = [LocalConfig.pymolTemplateName,'Ref_Pharmacophore',label]
      if owner.showProtein():
        show.append(LocalConfig.pymolProteinName)

      if hasattr(self,'volPs') and self.volPs:
        label = 'exvol'
        locs = [list(x) for x in self.volPs]
        colors = [SearchUtils.colors[-1]]*len(locs)
        owner.viewer3d.AddPharmacophore(locs,colors,label)
        show.append('exvol')
        
      for obj in show:
        owner.viewer3d.DisplayObject(obj)
    else:
      owner.viewer3d.ShowMol(self.mol,name=label,showOnly=False,
                             zoom=owner.centerView())
    if (LocalConfig.hBondsAlwaysEnabled or owner.showProtein()) \
       and owner.showHBonds():
      owner.displayHBonds(mol=label)
    if owner.showCollisions():
      owner.displayCollisions(mol=label)

    owner.viewer3d.SetDisplayUpdate(True)

  def compare(self,other,col,ascending):
    if not isinstance(other,EmbedItem):
      return 0
    else:
      if col==0:
        return cmp(self.keyVal,other.keyVal)
      elif col==1:
        return cmp(self.shapeScore,other.shapeScore)
      elif col==2:
        return cmp(self.rms,other.rms)
        
  def text(self,col):
    if col==0:
      return self.label
    elif col==1:
      return '%.2f'%self.shapeScore
    elif col==2:
      return '%.2f'%self.rms
    else:
      return '...'

  def getName(self):
    """ looks up the name of the molecule (which may have changed due to
    tree edits) and returns it
    """
    nm = self.parent().parent().getName()
    nm = '%s; mapping: %d; embedding: %d'%(nm,self.parent().keyVal,self.index)
    return nm
  
  def getMol(self):
    """ looks up the name of the molecule (which may have changed due to
    tree edits), sets the name property, and returns the molecule
    """
    self.mol.SetProp('_Name',self.getName())
    return self.mol

  def getLabel(self):
    return self.dispLabel


class MappingItem(MolAlignmentListItem):
  """ A class for storing a pharmacophore's mapping onto a molecule (collection of atom ids)
  """
  def __init__(self,mapping,mol,bounds,label,keyVal,*args):
    MolAlignmentListItem.__init__(self,label,keyVal,*args)
    self.mol = Chem.Mol(mol.ToBinary())
    if mol.HasProp("_Name"):
     self.mol.SetProp("_Name",mol.GetProp("_Name"))
    self.bounds = bounds
    self.activated = None
    self.mapping = mapping

    self.embeds = []
    self.setExpandable(True)
    self.numDesired=LocalConfig.numAlignmentsDesired
    self.maxToTry=LocalConfig.maxAlignmentsToTry
    self.contextMenu=None
    self.numGets=0
    self.initChemistry()

    self.highlights=[]
    for entry in self.mapping:
      self.highlights.append(entry)

    #print 'bounds for %s'%label
    #for row in self.bounds:
    #  print ' ',' '.join(['% 4.2f'%x for x in row])


  def initChemistry(self):
    self.mol.RemoveAllConformers()

    self.mol2d = Chem.MolFromMolBlock(Chem.MolToMolBlock(self.mol))
    self.mol2d.RemoveAllConformers()
    rdDepictor.Compute2DCoords(self.mol2d)
    Chem.Kekulize(self.mol2d)

    self.excludedVols = []
    owner = self.listView().owner
    if owner.excludedVols:
      for exVol in owner.excludedVols:
        origInfo = exVol.featInfo
        if len(origInfo)!=len(self.mapping): raise ValueError,'length mismatch'
        info = [None]*len(self.mapping)
        for i in range(len(self.mapping)):
          info[i] = (self.mapping[i],origInfo[i][1],origInfo[i][2])
        ev = ExcludedVolume.ExcludedVolume(info,exclusionDist=exVol.exclusionDist)
        ev.origPos = exVol.origPos
        self.excludedVols.append(ev)
      

  def generateEmbeds(self,numDesired=-1,randomSeed=-1):
    if numDesired<=0:
      numDesired = self.numDesired
    self.numGets+=1
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    owner = self.listView().owner
    dlg = qtUtils.ProgressDialogHolder('Generating Embeddings',self.numDesired)
    newItems = SearchUtils.GetMolEmbeddings(self.mol,self.bounds,self.mapping,
                                            owner.pcophore,
                                            maxToTry=self.maxToTry,
                                            numDesired=numDesired,
                                            excludedVolumes=self.excludedVols,
                                            randomSeed=randomSeed,
                                            useDirs=owner.useDirs,
                                            twoStageOptimization=True,
                                            progressCallback=dlg)
    if not newItems:
      QApplication.restoreOverrideCursor()
      qtUtils.warning(LocalConfig.resultsNoEmbeddingMessage)
      if not self.embeds:
        self.setEnabled(False)
        self.listView().tipper.add(self,
                                   LocalConfig.resultsNoEmbeddingMessage)
      return
    else:
      self.listView().tipper.remove(self)
      self.setEnabled(True)
    idxOffset=len(self.embeds)+1
    self.embeds.extend(newItems)
    self.embeds.sort()
    for i,embed in enumerate(newItems):
      e,embed,exVols = embed
      de = e-self.parent().energy()
      label='E=%.0f kcal/mol'%(de)
      if self.mol.HasProp("_Name"):
        embed.SetProp("_Name",self.mol.GetProp("_Name"))
      embed.SetProp('RelativeEnergy',str(de))
      itm = EmbedItem(self.mapping,embed,label,de,exVols,self,index=i+idxOffset)
      if de > self.listView().energyTol:
        itm.setEnabled(False)
    QApplication.restoreOverrideCursor()

  def getChildState(self):
    res = []
    child = self.firstChild()
    while child:
      embed = child.mol
      if embed.HasProp('_Name'):
        nm = embed.GetProp('_Name')
      else:
        nm = ''
      exVols = child.exVols
      de = child.keyVal
      rmsd = child.rms
      shapeScore = child.shapeScore
      label = child.label
      props = (nm,label,de,rmsd,shapeScore)
      res.append((embed,props,exVols))
      child = child.nextSibling()
    return res

  def restoreChildState(self,children):
    self.embeds = []
    for embed,props,exVols in children:
      nm,label,de,rmsd,shapeScore = props
      if nm: 
        embed.SetProp('_Name',nm)
      newItm = EmbedItem(self.mapping,embed,
                         label,de,exVols,
                         self,
                         index=len(self.embeds)+1,
                         preAligned=True)
      newItm.rms=rmsd
      newItm.shapeScore=shapeScore
      newItm.aligned=True
      newItm.mol.SetProp('AlignmentRMS',str(rmsd))
      newItm.mol.SetProp('AlignmentShapeScore',str(shapeScore))
      if de > self.listView().energyTol:
        newItm.setEnabled(False)
      self.embeds.append((de+self.parent().energy(),embed,exVols))

  def setOpen(self,state):
    if not self.embeds:
      self.generateEmbeds()
    QListViewItem.setOpen(self,state)


  def drawMol(self):
    owner = self.listView().owner
    owner.hitsMolCanvas.setMol(self.mol2d)
    for i,highlights in enumerate(self.highlights):
      color = SearchUtils.colors[i%len(SearchUtils.colors)]
      owner.hitsMolCanvas.highlightAtoms(highlights,
                                         highlightColor=color,
                                         highlightRadius=LocalConfig.featHighlightRadius,
                                         append=True)

  def getMoreAlignments(self,nToGet=-1):
    if nToGet<=0:
      nToGet = self.numDesired
    self.generateEmbeds(nToGet,randomSeed=23*(len(self.embeds)+self.numGets))
    

  def activate(self):
    if not self.isSelected():
      return
    QListViewItem.activate(self)
    self.drawMol()

  def contextMenuRequested(self,itm,pos,col):
    if self!=itm:
      return
    if not self.contextMenu:
      self.contextMenu = QPopupMenu(self.listView())
      self.contextMenu.insertItem('&More Alignments',self.getMoreAlignments)
    menu = self.contextMenu
    if menu:
      res = menu.exec_loop(pos)
    
class MolItem(MolAlignmentListItem):
  """ a class for storing a molecule in the alignment tree

  """
  def __init__(self,mol,bounds,label,keyVal,*args):
    MolAlignmentListItem.__init__(self,label,keyVal,*args)
    self._setMol(mol)
    self.bounds = bounds
    self.activated = None
    self.mappings = None
    self.setExpandable(True)
    self._energy = None
    self.contextMenu = None

  def _setMol(self,mol,is2d=False):
    self.mol = Chem.Mol(mol.ToBinary())
    if mol.HasProp("_Name"):
      self.mol.SetProp("_Name",mol.GetProp("_Name"))
    self.mol2d = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
    if not is2d:
      self.mol2d.RemoveAllConformers()
      rdDepictor.Compute2DCoords(self.mol2d)
    Chem.Kekulize(self.mol2d)

  def getChildState(self):
    res = []
    child = self.firstChild()
    while child:
      label = child.label
      keyVal = child.keyVal
      mapping = child.mapping
      bounds = child.bounds
      mol = child.mol
      if mol.HasProp('_Name'):
        nm = mol.GetProp('_Name')
      else:
        nm = ''
      childState = child.getChildState()
      props = (nm,label,mapping,bounds,keyVal)
      res.append((mol,props,childState))
      child = child.nextSibling()
    return res

  def restoreChildState(self,children):
    if children:
      nChildren=0
      self.mappings = []
      for mol,props,childState in children:
        nm,label,mapping,bounds,keyVal=props
        if nm:
          mol.SetProp('_Name',nm)
        nChildren += 1
        newItm = MappingItem(mapping,mol,bounds,label,nChildren,self)      
        newItm.restoreChildState(childState)
        self.mappings.append(mapping)
    else:
      self.mappings = None
    
  def contextMenuRequested(self,itm,pos,col):
    if self!=itm:
      return
    if not self.contextMenu:
      self.contextMenu = QPopupMenu(self.listView())
      self.contextMenu.insertItem('&Copy Molecule to Clipboad',self.copyToClipboard)
      self.contextMenu.insertItem('&Duplicate Molecule',self.duplicate)
      if hasCDX:
        self.contextMenu.insertItem('&Edit Molecule',self.edit)
    menu = self.contextMenu
    if menu:
      res = menu.exec_loop(pos)

  def copyToClipboard(self):
    try:
      smi = Chem.MolToSmiles(self.mol2d,kekuleSmiles=True)
    except:
      qtUtils.error('problems converting molecule to SMILES',exc_info=True)
    else:
      clip = QApplication.clipboard()
      clip.setText(smi)
    
  def duplicate(self):
    lv = self.listView()
    nm = 'Copy of '+self.mol.GetProp('_Name')
    newItm = MolItem(self.mol,self.bounds,nm,
                     0,lv)
    newItm.setVisible(1)
    lv.ensureItemVisible(newItm)

  def edit(self):
    if not hasCDX:
      qtUtils.logger.info('chemdraw not available')
      return
    try:
      chemdraw.CDXDisplay(Chem.MolToMolBlock(self.mol2d),
                          inFormat='chemical/x-mdl-molfile')
    except:
      qtUtils.warning('could not send molecule to chemdraw',exc_info=True)
      return
    ourWin = win32gui.GetActiveWindow()
    chemdraw.RaiseChemDraw()
    res = QMessageBox.information(None,"Make Edits",
                                  "Edit the molecule in Chemdraw.\nWhen you have finished, click the 'Ok' button.\nClick 'Cancel' to abort the edit.",
                                  QMessageBox.Ok,QMessageBox.Cancel)
    win32gui.SetForegroundWindow(ourWin)
    qtUtils.logger.debug('Result: %s'%str(res))
    if res==QMessageBox.Ok:
      try:
        mb = str(chemdraw.CDXGrab(outFormat='chemical/x-mdl-molfile'))
        mol = Chem.MolFromMolBlock(mb)
        bounds = rdDistGeom.GetMoleculeBoundsMatrix(mol)
      except:
        qtUtils.warning('could not grab or process molecule from chemdraw',exc_info=True)
        return
      else:
        origFragCount = len(Chem.GetMolFrags(self.mol2d))
        newFragCount = len(Chem.GetMolFrags(mol))
        if newFragCount>origFragCount:
          qtUtils.warning('molecule edit increased number of fragments from %d to %d.\nThis is dangerous.\nPlease check your molecule.'%(origFragCount,newFragCount))
      mol.SetProp('_Name',self.mol.GetProp('_Name'))
      self._setMol(mol)
      self.bounds = bounds
      self.startOver()

      self.listView().owner.hitsMolCanvas.setMol(self.mol2d)
      self.setEnabled(True)
    else:
      qtUtils.warning('Aborted')
  def startOver(self):
    if self.isOpen():
      self.setOpen(False)

    # remove our children:
    child = self.firstChild()
    while child:
      next = child.nextSibling()
      self.takeItem(child)
      child = next

    # clear our status variables:
    self.mappings = None
    self._energy=None
      
  def initChemistry(self,nAttempts=10,maxPasses=5):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    dlg = qtUtils.ProgressDialogHolder('Generating Baseline Conformer',nAttempts)
    self._energy = 1e8
    pkl = self.mol.ToBinary()
    for i in range(nAttempts):
      tmpMol = Chem.Mol(pkl)
      rdDistGeom.EmbedMolecule(tmpMol,randomSeed=23*(i+1))
      try:
        ff = ChemicalForceFields.UFFGetMoleculeForceField(tmpMol)
      except:
        pass
      else:
        needsMore = ff.Minimize()
        nPasses=0
        while needsMore and nPasses<maxPasses:
          needsMore = ff.Minimize()
          nPasses+=1
        if nPasses>=maxPasses:
          qtUtils.warning("too many passes required to optimize baseline conformer")
        e = ff.CalcEnergy()
        if e < self._energy:
          self._energy = e
          if self.mol.HasProp("_Name"):
            tmpMol.SetProp("_Name",self.mol.GetProp("_Name"))
          self.mol = tmpMol
        dlg(i+1)  
    QApplication.restoreOverrideCursor()
    #self.activated=True

  def energy(self):
    if self._energy is None:
      self.initChemistry()
    return self._energy

  def hasMapping(self):
    owner = self.listView().owner
    matches = SearchUtils.GetMolFeatMatches(self.mol,owner.featureFactory,owner.pcophore)
    if matches:
      failed,foo,foo,foo = EmbedLib.MatchPharmacophore(matches,self.bounds,
                                                       owner.pcophore,
                                                       useDownsampling=True,
                                                       use2DLimits=True,mol=self.mol,
                                                       useDirs=owner.useDirs)
    else:
      failed = True
    return len(matches),not failed

  def generateMappings(self):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    owner = self.listView().owner
    matches = SearchUtils.GetMolFeatMatches(self.mol,owner.featureFactory,owner.pcophore)
    sz = 1
    for row in matches: sz *= len(row)

    dlg = qtUtils.ProgressDialogHolder('Generating Matches',sz)
    #print Chem.MolToSmiles(self.mol)
    #print self.bounds.shape
    #print self.bounds[0]
    self.mappings  = EmbedLib.GetAllPharmacophoreMatches(matches,self.bounds,
                                                         owner.pcophore,
                                                         useDownsampling=1,
                                                         progressCallback=dlg,
                                                         use2DLimits=True,
                                                         mol=self.mol,
                                                         verbose=False)
    if not len(self.mappings):
      itm = QListViewItem(self)
      if matches:
        itm.setText(0,LocalConfig.resultsNoMappingMessage)
      else:
        itm.setText(0,LocalConfig.resultsMissingFeatsMessage)
      itm.setEnabled(False)
    else:
      for i in range(len(self.mappings)):
        mapping =[list(x.GetAtomIds()) for x in self.mappings[i]]
        label = 'Mapping-%d'%(i+1)
        itm = MappingItem(mapping,self.mol,self.bounds,label,i+1,self)
    QApplication.restoreOverrideCursor()

  def setup(self):
    QListViewItem.setup(self)
    
  def setOpen(self,state):
    self.activate()
    if self.mappings is None:
      self.generateMappings()
    QListViewItem.setOpen(self,state)
    
  def activate(self):
    if not self.isSelected():
      return
    QListViewItem.activate(self)
    self.listView().owner.hitsMolCanvas.setMol(self.mol2d)

  def okRename(self,col):
    MolAlignmentListItem.okRename(self,col)
    self.mol.SetProp('_Name',self.label)

class MolAlignmentList(QListView):
  """ the alignment list view itself

  """
  def __init__(self,owner,*args,**kwargs):
    self.energyTol = kwargs.get('energyTol',LocalConfig.embeddingEnergyTol)

    QListView.__init__(self,*args)
    self.addColumn('Molecules')
    self.addColumn('Shape Score')
    self.addColumn('Alignment RMSD')
    self.owner = owner
    self.setSorting(0)
    self.setShowSortIndicator(1)
    self.setDefaultRenameAction(QListView.Accept)
    #self.setResizeMode(QListView.LastColumn)
    self.setColumnWidthMode(0,QListView.Maximum)
    self.setRootIsDecorated(True)
    self.connect(self,
                 SIGNAL('currentChanged(QListViewItem *)'),
                 self.currentChanged)
    self.connect(self,
                 SIGNAL('contextMenuRequested(QListViewItem *,const QPoint &,int)'),
                 self.contextMenuRequested)
    self.connect(self,
                 SIGNAL('mouseButtonClicked(int, QListViewItem *,const QPoint &,int)'),
                 self.mouseButtonClicked)


    self.tipper = MolAlignmentListToolTip(self)
    
  def getAllInstancesOf(self,klass):
    res = []
    firstChild = self.firstChild()
    if firstChild:
      stack = [firstChild]
      while stack:
        head = stack.pop(0)
        if isinstance(head,klass):
          res.append(head)

        child = head.firstChild()
        if child:
          stack.insert(0,child)

        sib = head.nextSibling()
        if sib:
          stack.append(sib)
    return res
    
  def getCheckedEmbeds(self,attr='getMol'):
    leaves = self.getAllInstancesOf(EmbedItem)
    res = []
    for leaf in leaves:
      if hasattr(leaf,'isOn') and leaf.isOn():
        thing = getattr(leaf,attr)
        res.append(thing())
    return res

  def getMols(self,selectedOnly=False,attr='getMol'):
    res = []
    items = self.getAllInstancesOf(MolItem)
    for item in items:
      if not selectedOnly or item.isSelected():
        thing = getattr(item,attr)
        #try:
        thing = thing()
        #except TypeError:
        #  pass
        #print 'thing:',thing
        res.append(thing)
    return res
    

  def reinitializeChildren(self):
    child = self.firstChild()
    while child:
      child.startOver()
      #m,found = child.hasMapping()
      self.tipper.remove(child)
      #if not m or not found:
      #  child.setEnabled(False)
      #  if not m:
      #    self.tipper.add(child,LocalConfig.resultsMissingFeatsMessage)
      #  else:
      #    self.tipper.add(child,LocalConfig.resultsNoMappingMessage)
      #else:
      child.setEnabled(True)
      child = child.nextSibling()
    self.owner.hitsMolCanvas.reset()
    self.owner.cleanupViewer()
  def deleteItem(self,item):
    if not item:
      return
    parent = item.parent()
    if not parent: # it's top level
      parent = self
    nxt = item.nextSibling()
    if not nxt:
      nxt = item.itemAbove()
    parent.takeItem(item)
    if nxt:
      self.setSelected(nxt,True)
      nxt.activate()
    
    


  def mouseButtonClicked(self,button,item,pos,col):
    """

    We have to engage in this bit of hackery because Qt's QCheckListItems
    have the unfortunate property of setting their checked state whenever
    they are clicked on.  Here we check to see if we've clicked on a
    QCheckListItem and, if so, check to see if that click is in the box
    itself.

    """
    if button==1 and item and col==0 and isinstance(item,QCheckListItem):
      # check to see if we're clicking in the checkbox itself:
      mapped = self.mapFromGlobal(pos)
      x = mapped.x()
      d = item.depth()+1
      checkLeft = d*self.treeStepSize()+5
      checkRight = checkLeft+10
      if x<=checkRight and x>=checkLeft:
        item.setOn(not item.isOn())

  def currentChanged(self,item):
    if item:
        item.activate()
    else:
      self.owner.hitsMolCanvas.reset()
     
  def clear(self):
    QListView.clear(self)
    self.owner.hitsMolCanvas.reset()

  def contextMenuRequested(self,itm,pos,col):
    if hasattr(itm,'contextMenuRequested'):
      itm.contextMenuRequested(itm,pos,col)

  def setSorting(self,col,ascending=True):
    QListView.setSorting(self,col,ascending)
    # make sure that the currently selected objects remains
    # visible
    self.ensureItemVisible(self.selectedItem())

  def keyPressEvent(self,evt):
    handled = 0
    code = evt.key()
    item = None
    needRefresh = []
    if code == Qt.Key_Delete:
      item = self.selectedItem()
      if item is not None:
        self.deleteItem(item)
        handled=1
    if not handled:
      QListView.keyPressEvent(self,evt)
    else:
      evt.accept()

