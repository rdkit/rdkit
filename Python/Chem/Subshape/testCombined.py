import Chem
from Chem.PyMol import MolViewer
from Chem.Subshape import SubshapeBuilder,SubshapeObjects,SubshapeAligner
import cPickle,copy

m1 = Chem.MolFromMolFile('testData/square1.mol')
m2 = Chem.MolFromMolFile('testData/square2.mol')

b = SubshapeBuilder.SubshapeBuilder()
if 0:
  print 'm1:'
  s1 = b.GenerateSubshapeShape(m1)
  cPickle.dump(s1,file('testData/square1.shp.pkl','wb+'))
  print 'm2:'
  s2 = b.GenerateSubshapeShape(m2)
  cPickle.dump(s2,file('testData/square2.shp.pkl','wb+'))
else:
  s1 = cPickle.load(file('testData/square1.shp.pkl','rb'))
  s2 = cPickle.load(file('testData/square2.shp.pkl','rb'))

ns1 = copy.deepcopy(s1)
ns1.skelPts=None
b.GenerateSubshapeSkeleton(ns1)

v = MolViewer()
SubshapeObjects.DisplaySubshape(v,s1,'shape1')
SubshapeObjects.DisplaySubshape(v,ns1,'ns1')
#SubshapeObjects.DisplaySubshape(v,s2,'shape2')

a = SubshapeAligner.SubshapeAligner()
pruneStats={}
algs =a.GetSubshapeAlignments(m1,s1,m2,s2,b,pruneStats=pruneStats)
print len(algs)
print pruneStats

