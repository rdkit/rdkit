tests = [
  ("python", "UnitTestSDFToCSV.py", {}),
]

longTests = []

if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  doLong = 0
  if '-l' in sys.argv:
    doLong = 1
  failed, tests = TestRunner.RunScript('test_list.py', doLong, 1)
  sys.exit(len(failed))
