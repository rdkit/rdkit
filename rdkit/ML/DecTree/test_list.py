tests = [
  ("python", "UnitTestTreeUtils.py", {}),
  ("python", "UnitTestTree.py", {}),
  ("python", "UnitTestXVal.py", {}),
  ("python", "UnitTestID3.py", {}),
  ("python", "UnitTestPrune.py", {}),
  ("python", "UnitTestQuantTree.py", {}),
  ("python", "UnitTestSigTree.py", {}),
]

longTests = []
if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
