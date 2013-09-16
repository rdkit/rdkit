
tests=[
    ("python","testHelpers.py",{}),
    ]

longTests=[]

if __name__=='__main__':
  import sys
  import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
