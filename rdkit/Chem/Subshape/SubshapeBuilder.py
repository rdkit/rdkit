#
# Copyright (C) 2007 by Greg Landrum
#  All rights reserved
#

import copy
import pickle
import time

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem.Subshape import BuilderUtils, SubshapeObjects


class SubshapeCombineOperations(object):
  UNION = 0
  SUM = 1
  INTERSECT = 2


class SubshapeBuilder(object):
  gridDims = (20, 15, 10)
  gridSpacing = 0.5
  winRad = 3.0
  nbrCount = 7
  terminalPtRadScale = 0.75
  fraction = 0.25
  stepSize = 1.0
  featFactory = None

  def SampleSubshape(self, subshape1, newSpacing):
    ogrid = subshape1.grid
    rgrid = Geometry.UniformGrid3D(self.gridDims[0], self.gridDims[1], self.gridDims[2], newSpacing)
    for idx in range(rgrid.GetSize()):
      l = rgrid.GetGridPointLoc(idx)
      v = ogrid.GetValPoint(l)
      rgrid.SetVal(idx, v)

    res = SubshapeObjects.ShapeWithSkeleton()
    res.grid = rgrid
    return res

  def GenerateSubshapeShape(self, cmpd, confId=-1, addSkeleton=True, **kwargs):
    shape = SubshapeObjects.ShapeWithSkeleton()
    shape.grid = Geometry.UniformGrid3D(self.gridDims[0], self.gridDims[1], self.gridDims[2],
                                        self.gridSpacing)
    AllChem.EncodeShape(cmpd, shape.grid, ignoreHs=False, confId=confId)
    if addSkeleton:
      conf = cmpd.GetConformer(confId)
      self.GenerateSubshapeSkeleton(shape, conf, **kwargs)
    return shape

  def __call__(self, cmpd, **kwargs):
    return self.GenerateSubshapeShape(cmpd, **kwargs)

  def GenerateSubshapeSkeleton(self, shape, conf=None, terminalPtsOnly=False, skelFromConf=True):
    if conf and skelFromConf:
      pts = BuilderUtils.FindTerminalPtsFromConformer(conf, self.winRad, self.nbrCount)
    else:
      pts = BuilderUtils.FindTerminalPtsFromShape(shape, self.winRad, self.fraction)
    pts = BuilderUtils.ClusterTerminalPts(pts, self.winRad, self.terminalPtRadScale)
    BuilderUtils.ExpandTerminalPts(shape, pts, self.winRad)
    if len(pts) < 3:
      raise ValueError('only found %d terminals, need at least 3' % len(pts))

    if not terminalPtsOnly:
      pts = BuilderUtils.AppendSkeletonPoints(shape.grid, pts, self.winRad, self.stepSize)
    for pt in pts:
      BuilderUtils.CalculateDirectionsAtPoint(pt, shape.grid, self.winRad)
    if conf and self.featFactory:
      BuilderUtils.AssignMolFeatsToPoints(pts, conf.GetOwningMol(), self.featFactory, self.winRad)
    shape.skelPts = pts

  def CombineSubshapes(self, subshape1, subshape2, operation=SubshapeCombineOperations.UNION):
    cs = copy.deepcopy(subshape1)
    if operation == SubshapeCombineOperations.UNION:
      cs.grid |= subshape2.grid
    elif operation == SubshapeCombineOperations.SUM:
      cs.grid += subshape2.grid
    elif operation == SubshapeCombineOperations.INTERSECT:
      cs.grid &= subshape2.grid
    else:
      raise ValueError('bad combination operation')
    return cs


if __name__ == '__main__':  # pragma: nocover
  import pickle
  import tempfile

  from rdkit.Chem.PyMol import MolViewer

  # cmpd = Chem.MolFromSmiles('CCCc1cc(C(=O)O)ccc1')
  # cmpd = Chem.AddHs(cmpd)
  if 1:
    cmpd = Chem.MolFromSmiles('C1=CC=C1C#CC1=CC=C1')
    cmpd = Chem.AddHs(cmpd)
    AllChem.EmbedMolecule(cmpd)
    AllChem.UFFOptimizeMolecule(cmpd)
    AllChem.CanonicalizeMol(cmpd)
    # print(Chem.MolToMolBlock(cmpd), file=file('testmol.mol', 'w+'))
  else:
    cmpd = Chem.MolFromMolFile('testmol.mol')
  builder = SubshapeBuilder()
  if 1:
    shape = builder.GenerateSubshapeShape(cmpd)
  v = MolViewer()
  if 1:
    tmpFile = tempfile.NamedTemporaryFile(suffix='.grd', delete=False).name
    v.server.deleteAll()
    Geometry.WriteGridToFile(shape.grid, tmpFile)
    time.sleep(1)
    v.ShowMol(cmpd, name='testMol', showOnly=True)
    v.server.loadSurface(tmpFile, 'testGrid', '', 2.5)
  v.server.resetCGO('*')

  with open('subshape.pkl', 'w+') as f:
    pickle.dump(shape, f)
  for i, pt in enumerate(shape.skelPts):
    v.server.sphere(tuple(pt.location), .5, (1, 0, 1), 'Pt-%d' % i)
    if not hasattr(pt, 'shapeDirs'):
      continue
    momBeg = pt.location - pt.shapeDirs[0]
    momEnd = pt.location + pt.shapeDirs[0]
    v.server.cylinder(tuple(momBeg), tuple(momEnd), .1, (1, 0, 1), 'v-%d' % i)
