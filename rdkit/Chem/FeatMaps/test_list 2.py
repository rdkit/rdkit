tests = [
  ("python", "FeatMapPoint.py", {}),
  ("python", "UnitTestFeatMap.py", {}),
  ("python", "UnitTestFeatMapPoint.py", {}),
  ("python", "UnitTestFeatMapParser.py", {}),
  ("python", "UnitTestFeatMapUtils.py", {}),
]

longTests = []

if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
