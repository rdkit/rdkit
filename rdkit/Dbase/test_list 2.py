tests = [
  ("python", "UnitTestStorageUtils.py", {}),
  ("python", "UnitTestDbConnect.py", {}),
  ("python", "UnitTestDbInfo.py", {}),
  ("python", "UnitTestDbUtils.py", {}),
  ("python", "UnitTestDbResultSet.py", {}),
]

longTests = []

if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
