# $Id$
#
# Copyright (C) 2005 Rational Discovery LLC
#  All Rights Reserved
#
from qt import *
from qtGui import qtUtils
import Chem
from Chem.Pharm3D.ExcludedVolume import ExcludedVolume
from qtGui.Search3D import SearchUtils
from qtGui.Search3D.PharmacophoreDefWidgets import ExcludedVolumeWidget

# ------------------------------------------------------------------------
class ExcludedVolumeMixin:
  """

  We use this class to clean up the main 3D search widget code by collecting
  the excluded volume functionality in one place.

  This class cannot be instantiated and is intended solely to be mixed in to
  Search3DWidget.

  """
  def __init__(self,*args,**kwargs):
    raise TypeError,'Class cannot be instantiated'

  def initExcludeButtons(self):
    lout = QHBoxLayout()
    self.addExcludeButton=QPushButton(self.trUtf8('Add Excluded Atom(s)'),
                                      self.pcophorePage)
    lout.addWidget(self.addExcludeButton)
    self.connect(self.addExcludeButton,SIGNAL('clicked()'),
                 self.addExclusionClicked)
    QToolTip.add(self.addExcludeButton,self.trUtf8("Add the current atom picks as excluded volumes."))

    self.resetExcludeButton=QPushButton(self.trUtf8('Reset Excluded Atoms'),
                                        self.pcophorePage)
    lout.addWidget(self.resetExcludeButton)
    self.connect(self.resetExcludeButton,
                 SIGNAL('clicked()'),
                 self.resetExclusionClicked)
    QToolTip.add(self.resetExcludeButton,
                 self.trUtf8("Remove the current excluded volume settings."))
    self.pcophorePage.layout().addLayout(lout)

  def grabExcludedAtoms(self):
    """ returns a dictionary keyed by (label,idx) with coordinate values
      (label, idx) are the results from pymol

    """
    sels = self.viewer3d.GetSelectedAtoms('pkset')
    if not sels:
      sels = self.viewer3d.GetSelectedAtoms('pk1')
      if not sels:
        sels = self.viewer3d.GetSelectedAtoms()
    if not sels:
      qtUtils.warning("Please pick some atoms to be excluded")
      return None
    res = self.viewer3d.GetAtomCoords(sels)
    return res

  def addExcludedAtoms(self,excludedAtoms=None):
    if not excludedAtoms:
      excludedAtoms = self.grabExcludedAtoms()
    if not excludedAtoms:
      qtUtils.warning("No atoms provided")
      return
    self.excludedAtoms.update(excludedAtoms)

    featNames = [x.GetFamily() for x in self.activeFeats]

    if not hasattr(self,'exVolWidget') or not self.exVolWidget:
      self.exVolWidget = ExcludedVolumeWidget(self.pcophorePage,self,featNames)
      lout = self.pcophorePage.layout()
      lout.addWidget(self.exVolWidget)
      self.exVolWidget.show()

    for (label,idx),pos in excludedAtoms.iteritems():
      dists = []
      for feat in self.activeFeats:
        minD,maxD = SearchUtils.GetFeatToPointDistanceRange(self.mol,feat,pos)
        minD -= self.slop
        maxD += self.slop
        dists.append((minD,maxD))
      self.exVolWidget.addRow(label,idx,dists,pos)
      
  def resetExcludedAtoms(self):
    self.excludedAtoms = {}
    if hasattr(self,'exVolWidget'):
      self.exVolWidget.clearGrid()
      self.exVolWidget.setEnabled(False)
      self.exVolWidget.hide()
      self.exVolWidget.destroy()
    self.exVolWidget=None


  def buildExcludedVolumes(self):
    self.excludedVols = []
    if hasattr(self,'exVolWidget') and self.exVolWidget:
      for rowIdx,row in enumerate(self.exVolWidget.rows):
        featInfo = []
        for i,feat in enumerate(self.activeFeats):
          minV = float(str(row[i][0].text()))
          maxV = float(str(row[i][1].text()))
          ids = list(feat.GetAtomIds())
          featInfo.append((ids,minV,maxV))
        dist = float(str(row[len(self.activeFeats)].text()))
        ev = ExcludedVolume(featInfo,exclusionDist=dist)
        ev.feats = self.activeFeats[:]
        ev.origPos = self.exVolWidget.locs[rowIdx]
        self.excludedVols.append(ev)
    return self.excludedVols

  def enableExcludeButtons(self,state):
    if not hasattr(self,'addExcludeButton'):
      return
    if state:
      self.addExcludeButton.show()
      self.resetExcludeButton.show()
    else:
      self.addExcludeButton.hide()
      self.resetExcludeButton.hide()
