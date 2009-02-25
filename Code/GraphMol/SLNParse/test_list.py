import sys

tests=[
  ("testExecs/test.exe","",{}),
  ("python","test_list.py",{'dir':'Wrap'}),
  ]

if sys.platform != 'win32':
  pass



longTests=[
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
