tests = [
  ("python", "UnitTestFilter.py", {}),
  ("python", "UnitTestMLData.py", {}),
  ("python", "UnitTestQuantize.py", {}),
  ("python", "UnitTestStats.py", {}),
  ("python", "SplitData.py", {}),
  ("python", "DataUtils.py", {}),
]

longTests = []
if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
