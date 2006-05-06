from VLib import EnumNodes

o = EnumNodes.OutputNode(name='Out')
c = EnumNodes.CombineNode('C')
o.AddInput(c)
T = EnumNodes.SupplyNode([100,200,300],'Temp')
A = EnumNodes.AlloyNode('Alloy')
M1 = EnumNodes.SupplyNode(['Fe','Co','Ni'],'M1')
M2 = EnumNodes.SupplyNode(['Ru','Rh','Pd'],'M2')
A.AddInput(M1)
A.AddInput(M2)
S = EnumNodes.SupplyNode(['Support1','Support2','Support3'],'Support')
c.AddInput(A)
c.AddInput(S)
c.AddInput(T)


import cPickle
cPickle.dump(o,open('CatalEnum.pkl','wb+'))

def testf(*args):
  print args
f=EnumNodes.FilterNode(testf,'filter')
f.AddInput(A)

