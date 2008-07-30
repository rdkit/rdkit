# $Id$
#
# Copyright (C) 2005 Rational Discovery LLC
#  All Rights Reserved
#
from qt import *
from qtGui import qtUtils
import Chem
from Chem.Pharm3D.Pharmacophore import Pharmacophore
from qtGui.Search3D import SearchUtils,LocalConfig
from qtGui.Search3D.PharmacophoreDefWidgets import DistAssignmentWidget
import sets

# ------------------------------------------------------------------------
class PharmacophoreMixin:
  """

  We use this class to clean up the main 3D search widget code by collecting
  the pharmacophore functionality in one place.

  This class cannot be instantiated and is intended solely to be mixed in to
  Search3DWidget.

  """
  def __init__(self,*args,**kwargs):
    raise TypeError,'Class cannot be instantiated'
  def checkFDefs(self):
    " returns whether or not all feature types are defined "
    selections = []
    for row in self.featAssignmentWidget.rows:
      box = row[0]
      itm = box.currentItem()
      if itm<=0:
        label = str(row[2].text())
        label = label.replace(':','')
        qtUtils.error('No feature type selected for %s.'%label)
        return False
      selections.append(itm-1)
    self.selectedFeats=selections
    return True

  def grabDistances(self):
    if not self.checkFDefs():
      return
    if not self.expertModeEnabled() and \
       self.molAlignList.firstChild() is not None:
      confirm = QMessageBox.warning(self,
                              self.trUtf8("Discard Results?"),
                              self.trUtf8("Updating the pharmacophore will discard your current results. Continue?"),
                              1,2)
      if confirm != 1:
        return
    
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
    if self.distAssignmentWidget is None:
      self.distAssignmentWidget = DistAssignmentWidget(self.pcophorePage,self)
      lout = self.pcophorePage.layout()
      lout.insertWidget(5,self.distAssignmentWidget)
      space = QSpacerItem(10,10,QSizePolicy.Minimum,QSizePolicy.Expanding)
      lout.insertItem(6,space)
      self.distAssignmentWidget.show()
      self.distAssignmentWidget.setUpdatesEnabled(False)
    else:
      self.distAssignmentWidget.setUpdatesEnabled(False)
      self.distAssignmentWidget.reinit()

    nIndices = len(self.activeIndices)
    self.activeFeats = []
    assert nIndices==len(self.selectedFeats)
    for i,idx in enumerate(self.activeIndices):
      possibleFeats = self.featMap[idx-1]
      selected = possibleFeats[self.selectedFeats[i]]
      self.activeFeats.append(selected)

    # check for feature consistency:
    for i,feat in enumerate(self.activeFeats):
      atomsI = sets.Set(list(feat.GetAtomIds()))
      for j in range(i+1,len(self.activeFeats)):
        atomsJ = sets.Set(list(self.activeFeats[j].GetAtomIds()))
        if len(atomsI.intersection(atomsJ)):
          QApplication.restoreOverrideCursor()
          qtUtils.error('Features %d and %d have an atom in common.\nThis is not allowed.'%(i+1,j+1))
          return
    try:
      dists,dists2D = SearchUtils.GetFeatDistances(self.mol,self.activeFeats)
    except:
      QApplication.restoreOverrideCursor()
      qtUtils.error('A problem was encountered calculating feature--feature distances.',exc_info=True)
      return
    
    offset = 0
    for i in range(nIndices):
      featI = self.activeFeats[i]
      idxI = self.activeIndices[i]
      for j in range(i+1,nIndices):
        featJ = self.activeFeats[j]
        idxJ = self.activeIndices[j]
        currV = dists[offset]
        minV = currV - self.slop
        maxV = currV + self.slop

        max2D = dists2D[offset]
        max2D = int(max2D*(1+self.slop2D))
        
        self.distAssignmentWidget.addRow(featI,idxI,featJ,idxJ,minV,currV,maxV,max2D)
        offset += 1

    # this is a hacky work-around for a problem with qt and tables that
    # are more than 4 columns wide:
    hdr = self.distAssignmentWidget.verticalHeader()
    sz = 0
    for i in range(hdr.count()):
      sz += hdr.sectionSize(i)
    # and this is a hack on top of a hack: hardcode the header size:
    sz += 23
    self.distAssignmentWidget.setMaximumSize(QSize(1000,sz))
    self.distAssignmentWidget.setUpdatesEnabled(True)
    QApplication.restoreOverrideCursor()
      
    self.buildPharmacophore()
    if not self.expertModeEnabled():
      self.molAlignList.reinitializeChildren()
    else:
      self.enablePages(True,True,True)
      self.molAlignList.clear()

  def buildPharmacophore(self,
                         label=LocalConfig.refPcophoreName,
                         displayOn=LocalConfig.pymolTemplateName):
    if not hasattr(self,'activeFeats') or not self.activeFeats or \
       not hasattr(self,'distAssignmentWidget') or not self.distAssignmentWidget:
      qtUtils.error('Please define a pharmacophore before searching')
      return None
    self.pcophore = Pharmacophore(self.activeFeats)
    nFeats = len(self.activeFeats)
    offset = 0
    for i in range(nFeats):
      for j in range(i+1,nFeats):
        offset += 1
        row = self.distAssignmentWidget.rows[offset]
        minV = float(str(row[0].text()))
        maxV = float(str(row[1].text()))
        self.pcophore.setLowerBound(i,j,minV)
        self.pcophore.setUpperBound(i,j,maxV)
        max2D = row[2].value()
        self.pcophore.setUpperBound2D(i,j,max2D)
    self.displayPharmacophore(label=label,displayOn=displayOn)
    self.updateGeometry()
    self.pcophoreMolCanvas.updateGeometry()
    #print ' ------------- pcop draw'
    self.pcophoreMolCanvas.setMol(self.mol2d)
    self.displayPharmacophore2D()
    return self.pcophore

  def displayPharmacophore(self,
                           label=LocalConfig.refPcophoreName,
                           displayOn=''):
    if not hasattr(self,'pcophore') or not self.pcophore:
      qtUtils.error('no pharmacophore to display')
      return
    locs = [list(x.GetPos()) for x in self.pcophore.getFeatures()]
    colors = [SearchUtils.colors[x%(len(SearchUtils.colors))] for x in range(len(locs))]
    self.viewer3d.AddPharmacophore(locs,colors,label)
    if displayOn:
      self.viewer3d.Zoom(displayOn)

  def displayPharmacophore2D(self,highlightRad=LocalConfig.featHighlightRadius):
    if not hasattr(self,'activeFeats') or not self.activeFeats:
      qtUtils.error('no 2D pharmacophore to display (must have at least feature definitions)')
      return
    else:
      feats=self.activeFeats
      
    canvas = self.pcophoreMolCanvas
    canvas.clearHighlights()
    #print 'display2d'
    for i,feat in enumerate(feats):
      indices = feat.GetAtomIds()
      color = SearchUtils.colors[i]
      canvas.highlightAtoms(indices,highlightColor=color,
                            highlightRadius=highlightRad,append=True)
    canvas.update()




  
