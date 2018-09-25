tests = [
  ("python", "UnitTestBuildComposite.py", {}),
  ("python", "UnitTestScreenComposite.py", {}),
  ("python", "UnitTestAnalyzeComposite.py", {}),
]

for d in [
    'Cluster', 'Composite', 'Data', 'DecTree', 'Descriptors', 'InfoTheory', 'KNN', 'ModelPackage',
    'NaiveBayes', 'Neural', 'SLT', 'Scoring'
]:
  tests.append(('python', 'test_list.py', {'dir': d}))

longTests = []

if __name__ == '__main__':
  import sys
  from rdkit import TestRunner
  failed, tests = TestRunner.RunScript('test_list.py', 0, 1)
  sys.exit(len(failed))
