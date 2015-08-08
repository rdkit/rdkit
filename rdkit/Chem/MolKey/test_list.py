try:
  from rdkit.Chem.MolKey import MolKey
  tests=[
    ("python","InchiInfo.py",{}),
    ("python","MolKey.py",{}),
  ]
except ImportError:
  pass




longTests=[
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
