tests = [("python", "UnitTestCorrMatGen.py", {}), ]

longTests = [("python", "UnitTestBitRanker.py", {}), ]
if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
