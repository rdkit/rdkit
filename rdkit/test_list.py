
tests=[
  ("python","test_list.py",{'dir':'ML'}),
  ("python","test_list.py",{'dir':'Chem'}),
  ("python","test_list.py",{'dir':'DataStructs'}),
  ("python","test_list.py",{'dir':'Dbase'}),
  ("python","test_list.py",{'dir':'SimDivFilters'}),
  ("python","test_list.py",{'dir':'VLib'}),

  ]



longTests=[
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
