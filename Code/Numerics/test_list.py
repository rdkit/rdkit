
tests=[
  ("python","test_list.py",{'dir':'Optimizer'}),
  ("python","test_list.py",{'dir':'EigenSolvers'}),
  ("python","test_list.py",{'dir':'Alignment'}),
  ("testExecs/main.exe","",{}),

  ]


if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
