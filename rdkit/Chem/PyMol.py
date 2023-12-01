# $Id$
#
# Copyright (C) 2004-2012 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" uses pymol to interact with molecules

"""
import os
import sys
import tempfile

from rdkit import Chem

# Python3 compatibility
try:
  from xmlrpclib import Server
except ImportError:
  from xmlrpc.client import Server

_server = None


class MolViewer(object):

  def __init__(self, host=None, port=9123, force=0, **kwargs):
    global _server
    if not force and _server is not None:
      self.server = _server
    else:
      if not host:
        host = os.environ.get('PYMOL_RPCHOST', 'localhost')
      _server = None
      serv = Server('http://%s:%d' % (host, port))
      serv.ping()
      _server = serv
      self.server = serv
    self.InitializePyMol()

  def InitializePyMol(self):
    """ does some initializations to set up PyMol according to our
    tastes

    """
    self.server.do('set valence,1')
    self.server.do('set stick_rad,0.15')
    self.server.do('set mouse_selection_mode,0')
    self.server.do('set line_width,2')
    self.server.do('set selection_width,10')
    self.server.do('set auto_zoom,0')

  def DeleteAll(self):
    " blows out everything in the viewer "
    self.server.deleteAll()

  def DeleteAllExcept(self, excludes):
    " deletes everything except the items in the provided list of arguments "
    allNames = self.server.getNames('*', False)
    for nm in allNames:
      if nm not in excludes:
        self.server.deleteObject(nm)

  def LoadFile(self, filename, name, showOnly=False):
    """ calls pymol's "load" command on the given filename; the loaded object
    is assigned the name "name"
    """
    if showOnly:
      self.DeleteAll()
    id = self.server.loadFile(filename, name)
    return id

  def ShowMol(self, mol, name='molecule', showOnly=True, highlightFeatures=[], molB="", confId=-1,
              zoom=True, forcePDB=False, showSticks=False):
    """ special case for displaying a molecule or mol block """

    server = self.server
    if not zoom:
      self.server.do('view rdinterface,store')
    if showOnly:
      self.DeleteAll()

    if not forcePDB and mol.GetNumAtoms() < 999:
      if not molB:
        molB = Chem.MolToMolBlock(mol, confId=confId)
      mid = server.loadMolBlock(molB, name)
    else:
      if not molB:
        molB = Chem.MolToPDBBlock(mol, confId=confId)
      mid = server.loadPDB(molB, name)

    if highlightFeatures:
      nm = name + '-features'
      conf = mol.GetConformer(confId)
      for feat in highlightFeatures:
        pt = [0.0, 0.0, 0.0]
        for idx in feat:
          loc = conf.GetAtomPosition(idx)
          pt[0] += loc[0] / len(feat)
          pt[1] += loc[1] / len(feat)
          pt[2] += loc[2] / len(feat)
        server.sphere(pt, 0.2, (1, 1, 1), nm)
    if zoom:
      server.zoom('visible')
    else:
      self.server.do('view rdinterface,recall')
    if showSticks:  # show molecule in stick view
      self.server.do('show sticks, {}'.format(name))
    return mid

  def GetSelectedAtoms(self, whichSelection=None):
    " returns the selected atoms "
    if not whichSelection:
      sels = self.server.getNames('selections')
      if sels:
        whichSelection = sels[-1]
      else:
        whichSelection = None
    if whichSelection:
      items = self.server.index(whichSelection)
    else:
      items = []
    return items

  def SelectAtoms(self, itemId, atomIndices, selName='selection'):
    " selects a set of atoms "
    ids = '(id '
    ids += ','.join(['%d' % (x + 1) for x in atomIndices])
    ids += ')'
    cmd = 'select %s,%s and %s' % (selName, ids, itemId)
    self.server.do(cmd)

  def HighlightAtoms(self, indices, where, extraHighlight=False):
    " highlights a set of atoms "
    if extraHighlight:
      idxText = ','.join(['%s and (id %d)' % (where, x) for x in indices])
      self.server.do('edit %s' % idxText)
    else:
      idxText = ' or '.join(['id %d' % x for x in indices])
      self.server.do('select selection, %s and (%s)' % (where, idxText))

  def SetDisplayStyle(self, obj, style=''):
    " change the display style of the specified object "
    self.server.do('hide everything,%s' % (obj, ))
    if style:
      self.server.do('show %s,%s' % (style, obj))

  def SelectProteinNeighborhood(self, aroundObj, inObj, distance=5.0, name='neighborhood',
                                showSurface=False):
    """ selects the area of a protein around a specified object/selection name;
    optionally adds a surface to that """
    self.server.do('select %(name)s,byres (%(aroundObj)s around %(distance)f) and %(inObj)s' %
                   locals())

    if showSurface:
      self.server.do('show surface,%s' % name)
      self.server.do('disable %s' % name)

  def AddPharmacophore(self, locs, colors, label, sphereRad=0.5):
    " adds a set of spheres "
    self.server.do('view rdinterface,store')
    self.server.resetCGO(label)
    for i, loc in enumerate(locs):
      self.server.sphere(loc, sphereRad, colors[i], label, 1)
    self.server.do('enable %s' % label)
    self.server.do('view rdinterface,recall')

  def SetDisplayUpdate(self, val):
    if not val:
      self.server.do('set defer_update,1')
    else:
      self.server.do('set defer_update,0')

  def GetAtomCoords(self, sels):
    " returns the coordinates of the selected atoms "
    res = {}
    for label, idx in sels:
      coords = self.server.getAtomCoords('(%s and id %d)' % (label, idx))
      res[(label, idx)] = coords
    return res

  def HideAll(self):
    self.server.do('disable all')

  def HideObject(self, objName):
    self.server.do('disable %s' % objName)

  def DisplayObject(self, objName):
    self.server.do('enable %s' % objName)

  def Redraw(self):
    self.server.do('refresh')

  def Zoom(self, objName):
    self.server.zoom(objName)

  def DisplayHBonds(self, objName, molName, proteinName, molSelText='(%(molName)s)',
                    proteinSelText='(%(proteinName)s and not het)'):
    " toggles display of h bonds between the protein and a specified molecule "
    cmd = "delete %(objName)s;\n"
    cmd += "dist %(objName)s," + molSelText + "," + proteinSelText + ",mode=2;\n"
    cmd += "enable %(objName)s;"
    cmd = cmd % locals()

    self.server.do(cmd)

  def DisplayCollisions(self, objName, molName, proteinName, distCutoff=3.0, color='red',
                        molSelText='(%(molName)s)', proteinSelText='(%(proteinName)s and not het)'):
    " toggles display of collisions between the protein and a specified molecule "
    cmd = "delete %(objName)s;\n"
    cmd += "dist %(objName)s," + molSelText + "," + proteinSelText + ",%(distCutoff)f,mode=0;\n"
    cmd += """enable %(objName)s
    color %(color)s, %(objName)s"""
    cmd = cmd % locals()
    self.server.do(cmd)

  def SaveFile(self, filename):
    # PyMol will interpret the path to be relative to where it was started
    # from. Remedy that.
    if not filename:
      raise ValueError('empty filename')
    filename = os.path.abspath(filename)
    self.server.save(filename)

  def GetPNG(self, h=None, w=None, preDelay=0):
    try:
      import Image
    except ImportError:
      from PIL import Image
    import time
    if preDelay > 0:
      time.sleep(preDelay)
    fd = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
    fd.close()
    self.server.do('png %s' % fd.name)
    time.sleep(0.2)  # <- wait a short period so that PyMol can finish
    for i in range(10):
      try:
        img = Image.open(fd.name)
        break
      except IOError:
        time.sleep(0.1)
    try:
      os.unlink(fd.name)
    except (OSError, PermissionError):
      # happens sometimes on Windows. Not going to worry about this too deeply since
      # the files are in a temp dir anyway. This was github #936
      pass
    fd = None
    if h is not None or w is not None:
      sz = img.size
      if h is None:
        h = sz[1]
      if w is None:
        w = sz[0]
      if h < sz[1]:
        frac = float(h) / sz[1]
        w *= frac
        w = int(w)
        img = img.resize((w, h), True)
      elif w < sz[0]:
        frac = float(w) / sz[0]
        h *= frac
        h = int(h)
        img = img.resize((w, h), True)
    return img
