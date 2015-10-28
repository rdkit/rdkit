tests=[
  ("python","UnitTestcBitVect.py",{}),
  ("python","UnitTestcBitVect2.py",{}),
  ("python","UnitTestBitEnsemble.py",{}),
  ("python","UnitTestTopNContainer.py",{}),
  ("python","BitUtils.py",{}),
  ("python","VectCollection.py",{}),
  ("python","LazySignature.py",{}),
  ]
longTests=[
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
