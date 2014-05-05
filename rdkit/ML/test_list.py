
tests=[
  ("python","UnitTestBuildComposite.py",{}),
  ("python","UnitTestScreenComposite.py",{}),
  ("python","UnitTestAnalyzeComposite.py",{}),
  ]

for dir in ['Cluster','Composite','Data','DecTree','Descriptors','InfoTheory','KNN','ModelPackage','NaiveBayes','Neural','SLT']:
    tests.append(('python','test_list.py',{'dir':dir}))

longTests=[
  ]


if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
