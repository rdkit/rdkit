
tests=[
  ("python","SDSupply.py",{}),
  ("python","SmartsMolFilter.py",{}),
  ("python","SmartsRemover.py",{}),
  ("python","SmilesDupeFilter.py",{}),
  ("python","SmilesOutput.py",{}),
  ("python","SmilesSupply.py",{}),
  ]



longTests=[

  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
