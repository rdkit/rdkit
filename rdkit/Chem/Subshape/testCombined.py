import copy
import pickle

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.PyMol import MolViewer
from rdkit.Chem.Subshape import (SubshapeAligner, SubshapeBuilder,
                                 SubshapeObjects)

m1 = Chem.MolFromMolFile('test_data/square1.mol')
m2 = Chem.MolFromMolFile('test_data/square2.mol')

b = SubshapeBuilder.SubshapeBuilder()
b.gridDims = (10., 10., 5)
b.gridSpacing = 0.4
b.winRad = 2.0
if 1:
  print('m1:')
  s1 = b.GenerateSubshapeShape(m1)
  pickle.dump(s1, open('test_data/square1.shp.pkl', 'wb+'))
  print('m2:')
  s2 = b.GenerateSubshapeShape(m2)
  pickle.dump(s2, open('test_data/square2.shp.pkl', 'wb+'))
  ns1 = b.CombineSubshapes(s1, s2)
  b.GenerateSubshapeSkeleton(ns1)
  pickle.dump(ns1, open('test_data/combined.shp.pkl', 'wb+'))
else:
  s1 = pickle.load(open('test_data/square1.shp.pkl', 'rb'))
  s2 = pickle.load(open('test_data/square2.shp.pkl', 'rb'))
  #ns1 = pickle.load(file('test_data/combined.shp.pkl','rb'))
  ns1 = pickle.load(open('test_data/combined.shp.pkl', 'rb'))

v = MolViewer()
SubshapeObjects.DisplaySubshape(v, s1, 'shape1')
SubshapeObjects.DisplaySubshape(v, ns1, 'ns1')
#SubshapeObjects.DisplaySubshape(v,s2,'shape2')

a = SubshapeAligner.SubshapeAligner()
pruneStats = {}
algs = a.GetSubshapeAlignments(None, ns1, m1, s1, b, pruneStats=pruneStats)
print(len(algs))
print(pruneStats)

import os
import tempfile

from rdkit import Geometry

fName = tempfile.NamedTemporaryFile(suffix='.grd', delete=False).name
Geometry.WriteGridToFile(ns1.coarseGrid.grid, fName)
v.server.loadSurface(fName, 'coarse', '', 2.5)
os.unlink(fName)
fName = tempfile.NamedTemporaryFile(suffix='.grd', delete=False).name
Geometry.WriteGridToFile(ns1.medGrid.grid, fName)
v.server.loadSurface(fName, 'med', '', 2.5)
os.unlink(fName)
