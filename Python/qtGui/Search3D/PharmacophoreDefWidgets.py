# $Id$
#
# Copyright (C) 2005 Rational Discovery LLC
#  All Rights Reserved
#
from qt import *
from qttable import *
from qtGui.Search3D import SearchUtils
from qtGui import qtUtils

# ------------------------------------------------------------------------
class FeatAssignmentWidget(QWidget):
  def __init__(self,parent,searchWidget):
    QWidget.__init__(self,parent)
    self.searchWidget=searchWidget
    self.initGrid()

  def initGrid(self):
    self.gridLayout = QGridLayout(self,0,4)
    self.gridLayout.setColStretch(0,0)
    self.gridLayout.setColStretch(1,1)
    self.gridLayout.setColStretch(2,0)
    self.gridLayout.setColStretch(3,0)
    self.gridLayout.setSpacing(2)
    self.rows = []

  def removeRow(self,row):
    for entry in row:
      if isinstance(entry,QWidget):
        qtUtils.safeDestroyWidget(entry)
    
  def reinit(self):
    self.gridLayout.deleteAllItems()
    for row in self.rows:
      self.removeRow(row)
    self.rows=[]
    
  def addRow(self,id,feats):
    row = len(self.rows)
    txt = 'Atom %d'%id
    label = QLabel(txt+':',self)
    label.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
    self.gridLayout.addWidget(label,row,0)
    label.show()
    box = QComboBox(self)
    box.insertItem(' - Select Feature Type - ')
    for feat in feats:
      box.insertItem(feat.GetFamily())
    box.setEditable(False)
    box.show()
    boxcb = lambda x,y=feats:self.highlightFeat(x,y)
    self.connect(box,SIGNAL('highlighted(int)'),boxcb)
    self.gridLayout.addWidget(box,row,1)
    QToolTip.add(box,self.trUtf8("Set the feature type for %s"%txt))

    button = QPushButton("Highlight Atom",self)
    button.show()
    cb = lambda x=id:self.highlightAtom(x)
    self.connect(button,SIGNAL('clicked()'),cb)
    self.gridLayout.addWidget(button,row,2)
    QToolTip.add(button,self.trUtf8("Highlight atom %s in PyMol"%txt))

    button2 = QPushButton("Highlight Feature",self)
    button2.show()
    cb2 = lambda x=feats,y=box:self.highlightFeat(feats=x,box=y)
    self.connect(button2,SIGNAL('clicked()'),cb2)
    self.gridLayout.addWidget(button2,row,3)
    QToolTip.add(button,self.trUtf8("Highlight active feature %s in PyMol"%txt))

    self.gridLayout.setRowStretch(row,0)
    row = (box,boxcb,label,button,cb,button2,cb2,id)
    self.rows.append(row)
    return len(self.rows)

  def highlightAtom(self,id):
    self.searchWidget.highlightAtoms((id,),extraHighlight=True)

  def highlightFeat(self,which=-1,feats=None,box=None):
    if which<0 and box and hasattr(box,'currentItem'):
      which = box.currentItem()
    if which>0 and (feats and which-1<len(feats)):
      feat = feats[which-1]
      ids = [x+1 for x in feat.GetAtomIds()]
      self.searchWidget.highlightAtoms(ids,extraHighlight=False)
    else:
      self.searchWidget.clearAtomHighlights()
      self.searchWidget.update()

      
      
      
      

  def enableRow(self,idx,state):
    row = self.rows[idx]
    box = row[0]
    box.setEnabled(state)
    
  def enableRows(self,state,upTo=-1):
    while upTo>=0 and upTo<len(self.rows):
      self.enableRow(upTo,state)
      upTo-=1

  def isRowEnabled(self,row):
    box = row[0]
    return box.isEnabled()

  def getAtomIds(self):
    res = []
    for row in self.rows:
      res.append(row[-1])
    return tuple(res)

  def numRows(self):
    if self.rows:
      return len(self.rows)
    else:
      return 0

  def numEnabledRows(self):
    count = 0
    for row in self.rows:
      if self.isRowEnabled(row):
        count += 1
    return count

  def removeEnabledRows(self):
    newRows = []
    for row in self.rows:
      if self.isRowEnabled(row):
        self.removeRow(row)
      else:
        newRows.append(row)
    self.rows = newRows
    return len(newRows)



# ------------------------------------------------------------------------
class DistAssignmentWidget(QTable):
  def __init__(self,parent,searchWidget):
    QTable.__init__(self,parent)
    self.searchWidget=searchWidget
    self.colNames = ('Atom1','Atom2','Min Dist','Current Dist',
                     'Max Dist','Max 2D Dist','')
    self.stretchCols = (2,4)
    self.initGrid()
    #self.setSizePolicy(QSizePolicy(7,1,0,0,self.sizePolicy().hasHeightForWidth()))
  def initGrid(self):
    nCols = len(self.colNames)
    self.setNumCols(nCols)
    for i in range(nCols):
      self.setColumnStretchable(i,1)
    self.labelGrid()
    
  def labelGrid(self):
    row = []
    for i,thing in enumerate(self.colNames):
      self.horizontalHeader().setLabel(i,thing)
      row.append(thing)
    self.rows = [row]
    self.verticalHeader().hide()
    self.setLeftMargin(0)

  def reinit(self):
    for row in self.rows:
      for entry in row:
        if isinstance(entry,QWidget) and entry.parent()==self:
          qtUtils.safeDestroyWidget(entry)
    self.setNumRows(0)
    self.rows = []
    self.labelGrid()
    
  def addRow(self,feat1,id1,feat2,id2,minV,currV,maxV,max2D=10):
    row = len(self.rows)-1
    self.insertRows(row,1)
    txt1 = '%s-%d  '%(feat1.GetFamily(),id1)
    txt2 = '%s-%d  '%(feat2.GetFamily(),id2)
    l1 = QLabel(txt1,self)
    l2 = QLabel(txt2,self)
    l1.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
    l2.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
    self.setCellWidget(row,0,l1)
    self.setCellWidget(row,1,l2)
    l1.show()
    l2.show()

    e1 = QLineEdit('%.3f'%(minV),self)
    val1 = QDoubleValidator(0.0,1000,10,e1);
    e1.val=val1
    e1.setValidator(val1)
    l3 = QLabel('%.3f'%(currV),self)
    l3.setAlignment(Qt.AlignCenter|Qt.AlignVCenter)
    e2 = QLineEdit('%.3f'%(maxV),self)
    val2 = QDoubleValidator(0.0,1000,10,e2);
    e2.val=val2
    e2.setValidator(val2)
    self.setCellWidget(row,2,e1)
    self.setCellWidget(row,3,l3)
    self.setCellWidget(row,4,e2)
    e1.show()
    l3.show()
    e2.show()
    cb1 = lambda x,y=e1:self.textChanged(x,y)
    e1.cb = cb1
    cb2 = lambda x,y=e2:self.textChanged(x,y)
    e2.cb = cb2
    self.connect(e1,SIGNAL('textChanged(const QString &)'),cb1)
    self.connect(e2,SIGNAL('textChanged(const QString &)'),cb2)

    spin = QSpinBox(1,100,1,self)
    spin.setValue(max2D)
    spin.show()
    self.setCellWidget(row,5,spin)

    QToolTip.add(e1,
                 self.trUtf8("The minimum allowed distance between %s and %s."%(txt1,
                                                                                    txt2)))
    QToolTip.add(e2,
                 self.trUtf8("The maximum allowed distance between %s and %s."%(txt1,
                                                                                    txt2)))

    QToolTip.add(spin,
                 self.trUtf8("The maximum allowed 2D (topological) distance between %s and %s. This can help prevent highly folded conformations from matching."%(txt1,txt2)))

    b = QPushButton("Highlight",self)
    b.show()
    cb = lambda x=id1,y=id2:self.highlightClicked(x,y)
    self.connect(b,SIGNAL('clicked()'),cb)
    QToolTip.add(b,self.trUtf8("Highlight atoms %d and %d in PyMol"%(id1,id2)))
    self.setCellWidget(row,6,b)

    if not self.searchWidget.expertModeEnabled():
      e1.setEnabled(False)
      e2.setEnabled(False)
      spin.setEnabled(False)
    
    row = (e1,e2,spin,b,cb,l1,l2,l3)
    self.rows.append(row)
    self.adjustSize()
    self.setMinimumSize(self.minimumSizeHint())
    return len(self.rows)

  def enableRows(self,state,upTo=-1):
    while upTo>=0 and upTo<len(self.rows):
      row = self.rows[upTo]
      upTo-=1
      try:
        row[0].setEnabled(state)
        row[1].setEnabled(state)
        row[2].setEnabled(state)
      except AttributeError:
        pass
  def highlightClicked(self,id1,id2):
    self.searchWidget.highlightAtoms((id1,id2),extraHighlight=True)

  def textChanged(self,txt,which):
    if self.searchWidget.expertModeEnabled():
      self.searchWidget.enablePages(True,True,True)
      self.searchWidget.molAlignList.clear()

# ------------------------------------------------------------------------
class PharmacophoreTableItem(QTableItem):
  def __init__(self,owner,mol,match,bounds,*args):
    QTableItem.__init__(self,*args)
    self.owner = owner
    self.mol = mol
    self.match = match
    self.bounds = bounds

# ------------------------------------------------------------------------
class ExcludedVolumeWidget(QWidget):
  def __init__(self,parent,searchWidget,featNames):
    QWidget.__init__(self,parent)
    self.searchWidget=searchWidget
    self.featNames=featNames
    self.initGrid()

  def initGrid(self):
    nCols = 2*len(self.featNames)+3
    if not hasattr(self,'gridLayout') or not self.gridLayout:
      self.gridLayout = QGridLayout(self,0,nCols)

    for i in range(1,nCols-1):
      self.gridLayout.setColStretch(i,1)
    self.gridLayout.setSpacing(2)
    self.labelGrid()
    self.rows=[]
    self.locs=[]
    
  def labelGrid(self):
    row = []
    
    for i,featName in enumerate(self.featNames):
      col = 2*i+1
      label = QLabel(featName,self)
      label.setAlignment(Qt.AlignCenter|Qt.AlignVCenter)
      self.gridLayout.addMultiCellWidget(label,0,0,col,col+1)
      row.append(label)
      label.show()

    col = 2*len(self.featNames)+1
    label = QLabel('Radius',self)
    label.setAlignment(Qt.AlignCenter|Qt.AlignVCenter)
    self.gridLayout.setCellWidget(label,0,col)
    row.append(label)
    label.show()

    self.labelRows = row

  def clearGrid(self):
    self.gridLayout.deleteAllItems()
    for row in self.rows:
      for entry in row:
        if isinstance(entry,QWidget):
          qtUtils.safeDestroyWidget(entry)
    self.rows = []
    self.locs = []

  def reinit(self,featNames):
    self.clearGrid()
    self.featNames = featNames[:]
    self.initGrid()
    
  def addRow(self,where,id,currDists,pos,exclusionDist=3.0):
    rowIdx = len(self.rows)+1
    txt='Vol-%d'%(rowIdx)
    l1 = QLabel(txt,self)
    l1.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
    self.gridLayout.addWidget(l1,rowIdx,0)
    l1.show()

    row = []
    col = 1
    for i in range(len(self.featNames)):
      minV,maxV = currDists[i]
      e1 = QLineEdit('%.3f'%(minV),self)
      self.gridLayout.addWidget(e1,rowIdx,col)
      e2 = QLineEdit('%.3f'%(maxV),self)
      self.gridLayout.addWidget(e2,rowIdx,col+1)
      e1.show()
      e2.show()
      row.append( (e1,e2) )
      col+=2

      QToolTip.add(e1,
                   self.trUtf8("Set the minimum allowed distance between the %s feature and %s"%(self.featNames[i],txt)))
      QToolTip.add(e2,
                   self.trUtf8("Set the maximum allowed distance between the %s feature and %s"%(self.featNames[i],txt)))

    e1 = QLineEdit('%.3f'%(exclusionDist),self)
    self.gridLayout.addWidget(e1,rowIdx,col)
    e1.show()
    row.append(e1)
    QToolTip.add(e1,
                 self.trUtf8("Set the exclusion radius (closest allowed atom) for %s"%(txt)))

    col += 1
    button = QPushButton("Highlight",self)
    button.show()
    cb = lambda x=id,y=where:self.highlightClicked(x,y)
    self.connect(button,SIGNAL('clicked()'),cb)
    self.gridLayout.addWidget(button,rowIdx,col)
    row.append(button)
    row.append(cb)
    QToolTip.add(button,self.trUtf8("Highlight %s in PyMol"%txt))

    row.append(l1)

    self.rows.append(row)
    self.locs.append(pos)

  def highlightClicked(self,id,where):
    self.searchWidget.highlightAtoms((id,),where=where)
