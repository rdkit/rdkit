# $Id$
#
# Copyright (C) 2005-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" uses DSViewer to interact with molecules

"""
import os
import tempfile

from win32com.client import Dispatch

from rdkit import Chem

_nextDisplayId = 1


class Displayable(object):

  def __init__(self, doc, id=-1):
    global _nextDisplayId
    if id < 0:
      id = _nextDisplayId
      _nextDisplayId += 1

    self.doc = doc
    self.id = id
    self.visible = True
    self.children = []

  def Select(self, atoms=[], state=True, recurse=False):
    if state:
      selText = 'true'
    else:
      selText = 'false'
    if not atoms or atoms == '*':
      atomStr = '; atom "*"'
    else:
      # DSViewer has atom ids from 1, we do it from 0:
      atoms = ['id=%d' % (x) for x in atoms]
      atomStr = '; atom %s' % ','.join(atoms)

    cmd = 'SetProperty object RD_Visual=%d %s: select=%s' % (self.id, atomStr, selText)
    r = int(str(self.doc.DoCommand(cmd)))
    if not r and not atoms:
      # this handles an annoying case where if you try to select
      # a molecule by ID in DSViewer, you get nothing selected:
      atomStr = ''
      cmd = 'SetProperty object RD_Visual=%d %s: select=%s' % (self.id, atomStr, selText)
      r = int(str(self.doc.DoCommand(cmd)))
    #print 'sel cmd:',cmd
    #print 'result:', r

    # stupid DSViewer will select the bonds between pairs of highlighted atoms,
    # stop that nonsense right now:
    if r:
      cmd = 'SetProperty object RD_Visual=%d; bond index="*": select=off' % (self.id)
      self.doc.DoCommand(cmd)

    if recurse:
      for child in self.children:
        child.Select(atoms=atoms, state=state, recurse=True)
    return r

  def Hide(self, recurse=True):
    self.Select(state=True, recurse=True)
    self.doc.DoCommand('hide')
    self.Select(state=False, recurse=True)

  def Show(self, recurse=True):
    self.Select(state=True, recurse=True)
    self.doc.DoCommand('Show')
    self.Select(state=False, recurse=True)

  def ShowOnly(self, recurse=True):
    self.doc.DoCommand('HideAll')
    self.Select(state=True, recurse=True)
    self.doc.DoCommand('Show')
    self.Select(state=False, recurse=True)

  def __del__(self):
    self.doc.DoCommand('UnselectAll')
    count = self.Select(state=True, recurse=True)
    if count:
      self.doc.DoCommand('Delete')


class MolViewer(object):

  def __init__(self, force=0, title='Untitled', **kwargs):
    self.app = Dispatch('WebLabViewerPro.Application')
    self.app.Visible = 1
    if force or self.app.ActiveDocument is None:
      self.doc = self.app.New(title)
    else:
      self.doc = self.app.ActiveDocument

    self.displayables = {}

  def DeleteAll(self):
    self.doc.DoCommand('SelectAll')
    self.doc.DoCommand('Delete')
    self.displayables = {}

  def DeleteAllExcept(self, excludes):
    excludes = [x.lower() for x in excludes]
    allNames = self.displayables.keys()
    for nm in allNames:
      if nm not in excludes:
        del self.displayables[nm]

  def ShowMol(self, mol, name='molecule', showOnly=True, highlightFeatures=[], molB="", confId=-1,
              zoom=True):
    if showOnly:
      self.DeleteAll()
      obj = None
    else:
      obj = self.displayables.get(name.lower(), None)
      #if obj:
      #  obj.Select(state=True)
      #  self.doc.DoCommand('Delete')
      #  obj.Select(state=False)

    if not molB:
      molB = Chem.MolToMolBlock(mol, confId=confId)

    tmp = name + "\n" + molB[molB.index('\n') + 1:]
    molB = tmp

    if not obj:
      obj = Displayable(self.doc)
    if not hasattr(obj, '_molBlock') or obj._molBlock != molB:
      obj._molBlock = molB
      with tempfile.NamedTemporaryFile('w+', suffix='.mol') as tmpf:
        tmpf.write(molB)
      self.doc.DoCommand('PasteFrom %s' % tmp.name)
      self.doc.DoCommand('SetProperty molecule id=0 : RD_Visual=%d' % (obj.id))
      self.doc.DoCommand('SetProperty molecule id=0 : id=%d' % (obj.id))
      self.doc.DoCommand('SetProperty molecule id=0 : select=off')
    else:
      obj.Select(state=True)
      self.doc.DoCommand('Show')

    self.displayables[name.lower()] = obj

    if zoom:
      self.doc.DoCommand('Center')
      self.doc.DoCommand('FitView')

    return

  def LoadFile(self, filename, name, showOnly=False):
    if showOnly:
      self.DeleteAll()
    self.doc.DoCommand('PasteFrom %s' % filename)
    obj = Displayable(self.doc)
    self.doc.DoCommand('SetProperty molecule id=0 : id=%d' % (obj.id))
    self.doc.DoCommand('SetProperty molecule id=0 : select=off')
    count = self.doc.DoCommand('SetProperty AminoAcidChain id=0 : RD_Visual=%d' % (obj.id))
    if not count or int(count) <= 0:
      count = self.doc.DoCommand('SetProperty molecule id=0 : RD_Visual=%d' % (obj.id))
    self.displayables[name.lower()] = obj
    return obj

  def GetSelectedAtoms(self, whichSelection=''):
    #print 'WHICH',repr(whichSelection), whichSelection.lower() in self.displayables
    if not whichSelection:
      d = str(self.doc.DoCommand('GetPropertyValue atom select=true: id=?'))
      d2 = str(self.doc.DoCommand('GetPropertyValue atom select=true: molecule=?'))
      if d2:
        molIds = []
        tmpD = {}
        for id in d2.split(','):
          id = int(id.split('/')[1]) + 1
          if id in tmpD:
            molIds.append(tmpD[id])
          else:
            for k, v in self.displayables.iteritems():
              if id == v.id:
                tmpD[id] = k
                molIds.append(k)
      else:
        molIds = [''] * (d.count(',') + 1)
    elif whichSelection.lower() in self.displayables:
      whichSelection = whichSelection.lower()
      whichSelection = self.displayables[whichSelection].id
      d = str(
        self.doc.DoCommand('GetPropertyValue molecule RD_Visual=%d; atom select=true: id=?' %
                           whichSelection))
      molIds = [whichSelection] * (d.count(',') + 1)
    else:
      d = None
      molIds = None

    if d:
      splitD = d.split(',')
      #print 'splitD:',splitD
      #print 'molIds:',molIds
      try:
        res = []
        for i in range(len(splitD)):
          # DSViewer has atom ids from 1, we do it from 0:
          idx = int(splitD[i])
          res.append((molIds[i], idx))
      except Exception:
        import traceback
        traceback.print_exc()
        res = []
    else:
      res = []
    return res

  def HighlightAtoms(self, indices, where, extraHighlight=False):
    self.doc.DoCommand('UnSelectAll')
    self.SelectAtoms(where, indices)

  def SelectAtoms(self, itemId, atomIndices, selName='selection'):
    self.doc.DoCommand('UnSelectAll')
    self.doc.DoCommand('SetProperty atom id="*": select=off')
    o = self.displayables.get(itemId.lower(), None)
    #print 'O:',itemId,atomIndices
    if o:
      o.Select(atoms=atomIndices)

  def SetDisplayUpdate(self, val):
    if not val:
      self.doc.DoCommand('UpdateView off')
    else:
      self.doc.DoCommand('UpdateView on')

  def GetAtomCoords(self, sels):
    res = {}
    for label, idx in sels:
      whichSelection = label.lower()
      whichSelection = self.displayables[label].id
      # DSViewer has atom ids from 1, we do it from 0:
      idx += 1
      cmd = 'GetPropertyValue molecule RD_Visual=%d; atom id=%d: xyz=?' % (whichSelection, idx)
      coords = self.doc.DoCommand(cmd)
      coords = [float(x) for x in coords.split(' ')]
      res[(label, idx)] = coords
      #print 'grab:',label,idx,coords
    return res

  def AddPharmacophore(self, locs, colors, label, sphereRad=0.5):
    label = label.lower()
    self.SetDisplayUpdate(False)
    parent = Displayable(self.doc)
    for i, loc in enumerate(locs):
      color = colors[i]
      color = ' '.join([str(int(255 * x)) for x in color])
      obj = Displayable(self.doc)
      nm = 'sphere-%d' % obj.id
      self.doc.DoCommand('Sphere %s' % nm)
      self.doc.DoCommand('SetProperty Object name=%s : xyz=%f %f %f' % (nm, loc[0], loc[1], loc[2]))
      self.doc.DoCommand('SetProperty Object name=%s : radius=%f' % (nm, sphereRad))
      self.doc.DoCommand('SetProperty Object name=%s : color=%s' % (nm, color))
      self.doc.DoCommand('SetProperty Object name=%s : RD_Visual=%d' % (nm, parent.id))
      self.doc.DoCommand('SetProperty Object name=%s : id=%d' % (nm, parent.id))
      #parent.children.append(obj)
    self.displayables[label] = parent
    self.SetDisplayUpdate(True)

  def SetDisplayStyle(self, obj, style=''):
    self.doc.DoCommand('UnSelectAll')
    obj = obj.lower()
    o = self.displayables.get(obj, None)
    if o:
      o.Select(state=True)
      if style == 'sticks':
        self.doc.DoCommand('DisplayStyle Atom Stick')
      elif style == 'lines':
        self.doc.DoCommand('DisplayStyle Atom Line')
      elif style == '':
        self.doc.DoCommand('DisplayStyle Atom Off')
      o.Select(state=False)

  def HideAll(self):
    self.doc.DoCommand('HideAll')

  def HideObject(self, objName):
    self.doc.DoCommand('UnSelectAll')
    objName = objName.lower()
    o = self.displayables.get(objName, None)
    if o:
      o.Hide()

  def DisplayObject(self, objName):
    self.doc.DoCommand('UnSelectAll')
    objName = objName.lower()
    o = self.displayables.get(objName, None)
    if o:
      o.Show()

  def Zoom(self, objName):
    self.doc.DoCommand('UnSelectAll')
    objName = objName.lower()
    o = self.displayables.get(objName, None)
    if o:
      r = o.Select(state=True)
      self.doc.DoCommand('Center')
      self.doc.DoCommand('FitView')
      o.Select(state=False)

  def SelectProteinNeighborhood(self, aroundObj, inObj, distance=5.0, name='neighborhood',
                                showSurface=False):
    """ FIX: the surface display stuff here is all screwed up due to
    differences between the way PyMol and DSViewer handle surfaces.
    In PyMol they are essentially a display mode for the protein, so
    they don't need to be managed separately.
    In DSViewer, on the other hand, the surface is attached to the
    protein, but it needs to be hidden or shown on its own.  I haven't
    figured out how to do that yet.
    """
    self.doc.DoCommand('UnSelectAll')
    o = self.displayables.get(aroundObj.lower(), None)
    p = self.displayables.get(inObj.lower(), None)
    if o and p:
      self.SetDisplayUpdate(False)
      p.Show()
      self.doc.DoCommand('UnSelectAll')
      tmp = self.doc.DoCommand('SetProperty object RD_Visual=%d;object id="*":select=on' % o.id)
      tmp = self.doc.DoCommand('SelectByRadius inside %f atom' % distance)
      # that selects all atoms in the radius, now we need to make sure
      #  only atoms in _inObj_ are selected:
      for obj in self.displayables.values():
        if obj.id != p.id:
          self.doc.DoCommand('SetProperty object RD_Visual=%d;object id="*":select=off' % obj.id)

      # ----
      # now get all the residue names for the selected atoms:
      rs = self.doc.DoCommand('GetPropertyValue atom select=true: parent=?')
      if rs:
        rs = rs.split(',')
        residues = {}
        for r in rs:
          residues[r] = 1

        # and select each atom in those residues:
        parents = ','.join(['parent="%s"' % x for x in residues.keys()])
        cmd = 'SetProperty atom %s: select=on' % parents
        tmp = self.doc.DoCommand(cmd)
      if showSurface:
        # create the surface:
        self.doc.DoCommand('Surface')
        obj = Displayable(self.doc)
        self.displayables[name] = obj
        self.doc.DoCommand('SetProperty surface id="*":RD_Visual=%d' % obj.id)

        self.doc.DoCommand('UnSelectAll')

      self.SetDisplayUpdate(True)

  def Redraw(self):
    self.SetDisplayUpdate(True)


if __name__ == '__main__':
  from rdkit import Chem
  from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

  m = Chem.MolFromSmiles('c1cccc2c1cccc2')
  rdDistGeom.EmbedMolecule(m)
  rdForceFieldHelpers.UFFOptimizeMolecule(m)

  s = MolViewer()
  s.ShowMol(m)
