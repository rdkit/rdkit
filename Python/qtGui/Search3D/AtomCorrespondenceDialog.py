# $Id$
#
#  Copyright (C) 2006 Greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
from qt import *
from qttable import *
from forms.AtomCorrespondenceDialog import AtomCorrespondenceDialog as _Form
from qtGui.Search3D import LocalConfig,SearchUtils
from qtGui import qtUtils
from copy import deepcopy
import Chem
from Chem import ChemicalForceFields,rdShapeHelpers
from Numeric import *

class LocalComboTableItem(QComboTableItem):
  " we just use these to get a bit more space in tables "
  margin=(0,4)
  def sizeHint(self):
    sz = QComboTableItem.sizeHint(self)
    res = QSize(sz.width()+self.margin[0],sz.height()+self.margin[1])
    return res

class LocalTableItem(QTableItem):
  def __init__(self,tbl,text):
    QTableItem.__init__(self,tbl,QTableItem.Never,text)
    self.setText(text)
  def key(self):
    v = self.getInt()
    return 'X'*v

  def setText(self,text):
    if text=='0':
      text='None'
    QTableItem.setText(self,text)

  def setInt(self,v):
    if v<=0:
      text='None'
    else:
      text = str(v)
    QTableItem.setText(self,text)

  def getInt(self):
    try:
      res = int(str(self.text()))
    except:
      res = 0
    return res

class LocalTable(QTable):
  def sortColumn(self,col,ascending,wholeRows):
    QTable.sortColumn(self,col,ascending,True)

class AtomCorrespondenceDialog(_Form):
  def __init__(self,embedItem,refMol,parent=0,name='',modal=False,flags=0):
    _Form.__init__(self,parent,name,modal,flags)
    self.embedItem=embedItem
    self.setProbeMol(embedItem.getMol())
    self.refMol=refMol

    self.idTable = LocalTable(self,'idTable')
    self.idTable.setNumRows(0)
    self.idTable.setNumCols(3)
    self.idTable.horizontalHeader().setLabel(0,"Alignment Id")
    self.idTable.horizontalHeader().setLabel(1,"Reference Id")
    self.idTable.horizontalHeader().setLabel(2,"Button")
    self.idTable.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.idTable.sizePolicy().hasHeightForWidth()))
    self.idTable.setMinimumSize(QSize(0,0))
    self.idTable.setSelectionMode(QTable.Single)
    self.idTable.setColumnReadOnly(0,True)
    self.idTable.setColumnReadOnly(1,True)
    self.idTable.setColumnStretchable(2,True)
    self.idTable.horizontalHeader().setLabel(2,'')
    self.idTable.setSorting(True)
    self.layout().insertWidget(0,self.idTable)

    
    self.initCorrespondence()
    self.connect(self.idTable,SIGNAL('currentChanged(int,int)'),self.currentChanged)
    self.connect(self.idTable,SIGNAL('valueChanged(int,int)'),self.valueChanged)
    self.setIcon(QPixmap(qtUtils.logoImageData))
    
  def setProbeMol(self,mol):
    self.probeMol=deepcopy(mol)
    if mol.HasProp('_Name'):
      self.probeMol.SetProp('_Name',mol.GetProp('_Name'))

  def show(self):
    self.drawCurrentState()
    _Form.show(self)

  def getAlignmentProps(self):
    ff = ChemicalForceFields.UFFGetMoleculeForceField(self.probeMol)
    e = ff.CalcEnergy()
    de = e-self.embedItem.parent().parent().energy()

    owner = self.embedItem.listView().owner
    probeFeats = SearchUtils.GetFeaturesForAtomMatch(self.probeMol,self.embedItem.mapping,
                                                     owner.pcophore,owner.featureFactory)
    refFeats = owner.activeFeats
    rmsd = SearchUtils.GetAlignmentRMSD(self.refMol,refFeats,self.probeMol,probeFeats,
                                        useDirs=owner.useDirs)

    tani = rdShapeHelpers.ShapeTanimotoDist(self.probeMol,
                                            self.refMol,
                                            gridSpacing=0.5)
    shapeScore = 1.0-tani

    return de,shapeScore,rmsd
  
  def accept(self):
    from qtGui.Search3D.MolAlignmentList import EmbedItem
    index = self.embedItem.parent().childCount()+1

    de,shapeScore,rmsd = self.getAlignmentProps()
    label='E=%.0f kcal/mol'%(de)

    exVols = []
    newItm = EmbedItem(self.embedItem.mapping,self.probeMol,
                       label,de,exVols,
                       self.embedItem.parent(),preAligned=True)
    newItm.rms=rmsd
    newItm.shapeScore=shapeScore
    newItm.aligned=True
    newItm.mol.SetProp('AlignmentRMS',str(rmsd))
    newItm.mol.SetProp('AlignmentShapeScore',str(shapeScore))
    newItm.listView().ensureItemVisible(newItm)
    newItm.listView().setSelected(newItm,True)

    self.hide()
  def reject(self):
    self.hide()

  def initCorrespondence(self,noDupes=LocalConfig.noDupeCorrespondences,
                         maxDist=LocalConfig.maxAtomAtomCorrespondenceDistance):
    if not self.probeMol or not self.refMol:
      raise ValueError,'missing molecule'
    nAtoms = self.probeMol.GetNumAtoms()
    self.idTable.setNumRows(nAtoms)
    corresp = SearchUtils.GuessMolMolCorrespondence(self.refMol,self.probeMol,
                                                    noDupes=noDupes,maxDist=maxDist)
    for i in range(nAtoms):
      self.idTable.setItem(i,0,LocalTableItem(self.idTable,str(i+1)))
      sel = corresp.get(i,-1)+1
      self.idTable.setItem(i,1,LocalTableItem(self.idTable,str(sel)))

      cb = lambda atomId=i+1:self.clearCorrespondence(atomId)
      button = QPushButton('Clear Current',self.idTable)
      QToolTip.add(button,'Clears the correspondence for this atom.')
      self.connect(button,SIGNAL('clicked()'),cb)
      button.cb = cb
      self.idTable.setCellWidget(i,2,button)

    self.idTable.adjustColumn(0)
    self.idTable.adjustColumn(1)
    # add a bit of extra space to accomodate bold fonts:
    self.idTable.setColumnWidth(0,self.idTable.columnWidth(0)+5)
    self.idTable.setColumnWidth(1,self.idTable.columnWidth(1)+5)

  def findAtomRow(self,atomId):
    row = 0
    while row < self.idTable.numRows() and \
          self.idTable.item(row,0).getInt() != atomId:
      row+=1
    if row==self.idTable.numRows():
      raise ValueError,'bad atom id'
    return row

  def getPicks(self):
    owner = self.embedItem.listView().owner
    viewer = owner.viewer3d
    sel = viewer.GetSelectedAtoms('pkset')
    if not sel:
      sel = viewer.GetSelectedAtoms('pk1')
    if not sel:
      sel = viewer.GetSelectedAtoms()
    return sel

  def setCorrespondence(self,row,idx,validate=True):
    self.idTable.item(row,1).setInt(idx)
    while self.idTable.numSelections():
      self.idTable.removeSelection(0)
    sel = QTableSelection()
    sel.init(row,0)
    sel.expandTo(row,0)
    self.idTable.ensureVisible(row,0)
    self.idTable.addSelection(sel)
    if validate:
      self.validateCorrespondence()

  def assignDoubleSelection(self):
    sel = self.getPicks()
    
    refSel=[]
    probeSel=[]
    for item,idx in sel:
      if item==LocalConfig.pymolTemplateName:
        refSel.append(idx)
      elif item==LocalConfig.refinedAlignmentName:
        probeSel.append(idx)
    if not len(refSel):
      qtUtils.error('You must select an atom in the reference molecule.')
      return
    if len(refSel)!=1:
      qtUtils.error('You must select a single atom in the reference molecule.')
      return
    if not len(probeSel):
      qtUtils.error('You must select one atom in the probe alignment.')
      return
    if len(probeSel)!=1:
      qtUtils.error('You may select only one atom in the probe alignment.')
      return

    row = self.findAtomRow(probeSel[0])
    idx = refSel[0]
    self.setCorrespondence(row,idx)
    self.showCorrespondence()

  def clearCorrespondence(self,atomId):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    row = self.findAtomRow(atomId)
    self.idTable.item(row,1).setInt(0)
    self.idTable.updateCell(row,1)
    self.showCorrespondence()
    QApplication.restoreOverrideCursor()

  def drawCurrentState(self):
    owner = self.embedItem.listView().owner
    viewer = owner.viewer3d
    label=LocalConfig.refinedAlignmentName
    viewer.HideAll()
    viewer.ShowMol(self.probeMol,name=label,showOnly=False,
                   zoom=owner.centerView())
    show = [LocalConfig.pymolTemplateName,'Ref_Pharmacophore',label]
    if owner.showProtein():
      show.append(LocalConfig.pymolProteinName)
    for obj in show:
      owner.viewer3d.DisplayObject(obj)
    self.showCorrespondence()
    
  def getCurrentCorrespondence(self):
    corresp = {}
    for i in range(self.idTable.numRows()):
      row = self.idTable.item(i,0).getInt()-1
      sel = self.idTable.item(i,1).getInt()
      if sel>0:
        corresp[row]=sel-1
    return corresp

  def refineClicked(self):
    if not self.probeMol or not self.refMol or not self.embedItem:
      raise ValueError,'missing value'
    corresp = self.getCurrentCorrespondence()
    if not corresp:
      return
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    try:
      perf=SearchUtils.ImproveMolMolAlignment(self.refMol,self.probeMol,corresp)
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.warning('Problems encountered while improving alignment',
                      exc_info=True)
    else:
      de,shapeScore,rmsd = self.getAlignmentProps()
      self.energyBox.setText('%.0f'%de)
      self.shapeScoreBox.setText('%.2f'%shapeScore)
      self.rmsdBox.setText('%.2f'%rmsd)
      QApplication.restoreOverrideCursor()
      self.drawCurrentState()


  def validateCorrespondence(self):
    refAtomsSeen={}
    for i in range(self.idTable.numRows()):
      row = self.idTable.item(i,0).getInt()-1
      sel = self.idTable.item(i,1).getInt()
      if sel>0:
        refIdx = sel-1
        l = refAtomsSeen.get(refIdx,[])
        l.append(row)
        refAtomsSeen[refIdx]=l

    for refIdx,rowL in refAtomsSeen.iteritems():
      if len(rowL)>1:
        viewer=self.embedItem.listView().owner.viewer3d
        idTxt = ','.join([str(x+1) for x in rowL])
        cmd = """
        unpick
        delete offending_atoms
        select offending_atoms,(%s and id %d) or (%s and id %s)
        """%(LocalConfig.pymolTemplateName,refIdx+1,
             LocalConfig.refinedAlignmentName,idTxt)
        viewer.server.do(cmd)
        qtUtils.error('More than one probe atom corresponds to a single reference atom.\nPlease remove all but one of these correspondences.\nThe offending atoms (%s) are highlighted (selected) in PyMol.'%(idTxt))
        return False
    return True

  def valueChanged(self,row,col):
    if col==1:
      self.validateCorrespondence()
    
  def currentChanged(self,row,col):
    owner = self.embedItem.listView().owner
    viewer = owner.viewer3d
    # Note: we can constantly redraw here and make things
    # a bit more robust, but that slows things down enormously
    #self.drawCurrentState()

    atomId = self.idTable.item(row,0).getInt()

    # FIX: This is currently PyMol specific
    cmd = 'edit %s and (id %d)'%(LocalConfig.refinedAlignmentName,
                                 atomId)
                                                 
    itm = self.idTable.item(row,1).getInt()
    if itm>=0:
      cmd += ', %s and (id %d)'%(LocalConfig.pymolTemplateName,itm)
    viewer.server.do(cmd)

  def makeCorrespondenceClicked(self):
    self.assignDoubleSelection()

  def showCorrespondence(self):
    viewer = self.embedItem.listView().owner.viewer3d
    correspName=LocalConfig.pymolCorrespondenceName
    viewer.server.do('delete %(correspName)s'%locals())
    if self.showCorrespondenceBox.isChecked():
      refName=LocalConfig.pymolTemplateName
      probeName=LocalConfig.refinedAlignmentName
      for i in range(self.idTable.numRows()):
        probeIdx = self.idTable.item(i,0).getInt()
        sel = self.idTable.item(i,1).getInt()
        if sel:
          viewer.server.do('dist %(correspName)s,\
             (%(probeName)s and id %(probeIdx)d),\
             (%(refName)s and id %(sel)d),mode=0'%locals())

  def showCorrespondenceChanged(self,state):
    self.showCorrespondence()

  def clearCorrespondenceClicked(self):
    confirm = QMessageBox.warning(self,
                        self.trUtf8("Clear Correspondence?"),
                        self.trUtf8("Are you sure that you wish to clear the current atom-atom correspondence?"),
                        1,2)
    if confirm != 1:
      return
    
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    for row in range(self.idTable.numRows()):
      self.idTable.item(row,1).setInt(0)
      self.idTable.updateCell(row,1)
    self.showCorrespondence()
    QApplication.restoreOverrideCursor()

      
