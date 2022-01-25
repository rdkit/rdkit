import sys
from rdkit import TestRunner

tests = [
  ("python", "Utils.py", {}),
  ("python", "Torsions.py", {}),
  ("python", "Pairs.py", {}),
  ("python", "UnitTestDescriptors.py", {}),
]

longTests = []

if __name__ == '__main__':
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
