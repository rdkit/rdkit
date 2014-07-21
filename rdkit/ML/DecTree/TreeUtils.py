#
#  Copyright (C) 2001-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" Utilities for working with trees

"""

def CollectLabelLevels(tree,levels,level=0,maxDepth=1e8):
  if level < maxDepth:
    if not tree.GetTerminal():
      l = tree.GetLabel()
      currLevel = levels.get(l,1e8)
      if level < currLevel:
        levels[l] = level
      for child in tree.GetChildren():
        CollectLabelLevels(child,levels,level+1,maxDepth)
  return levels
  
def CollectDescriptorNames(tree,names,level=0,maxDepth=1e8):
  if level < maxDepth:
    if not tree.GetTerminal():
      names[tree.GetLabel()] = tree.GetName()
      for child in tree.GetChildren():
        CollectDescriptorNames(child,names,level+1,maxDepth)
  return names    
    
#------------------------------------
#
#  doctest boilerplate
#
_test1="""
>>> from DecTree import DecTreeNode as Node
>>> t1 = Node(None,'d1',1)
>>> t2 = Node(None,'d2',2)
>>> t1.AddChildNode(t2)
>>> t2 = Node(None,'d3',3)
>>> t1.AddChildNode(t2)
>>> t3 = Node(None,'d4',4)
>>> t2.AddChildNode(t3)
>>> t3 = Node(None,'d2',2)
>>> t2.AddChildNode(t3)
>>> r = CollectLabelLevels(t1,{})
>>> r[2]
1
>>> r[1]
0
>>> r[3]
1
>>> r[4]
2
>>> r = CollectLabelLevels(t1,{},0,2)
>>> r[2]
1
>>> r[1]
0
>>> r[3]
1
>>> 4 in r
0

Check that we can handle subtrees:
>>> r = CollectLabelLevels(t1,{},1,2)
>>> r[1]
1
>>> 2 in r
0
>>> 3 in r
0
>>> 4 in r
0

>>> names = CollectDescriptorNames(t1,{})
>>> names[1]
'd1'
>>> names[2]
'd2'
>>> names[3]
'd3'
>>> names[4]
'd4'

>>> names = CollectDescriptorNames(t1,{},0,2)
>>> names[1]
'd1'
>>> names[2]
'd2'
>>> names[3]
'd3'
>>> 4 in names
0

>>> names = CollectDescriptorNames(t1,{},1,2)
>>> names[1]
'd1'
>>> 2 in names
0
>>> 3 in names
0
>>> 4 in names
0


"""

__test__={'_test1':_test1}
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)

