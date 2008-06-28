
tests=[
  ("python","test_list.py",{'dir':'RDGeneral'}),
  ("python","test_list.py",{'dir':'DataStructs'}),
  ("python","test_list.py",{'dir':'DataStructs/Wrap'}),
  ("python","test_list.py",{'dir':'Query'}),
  ("python","test_list.py",{'dir':'DataManip/MetricMatrixCalc'}),
  ("python","test_list.py",{'dir':'DataManip/MetricMatrixCalc/Wrap'}),
  ("python","test_list.py",{'dir':'ML/InfoTheory/Wrap'}),
  ("python","test_list.py",{'dir':'Numerics'}),
  ("python","test_list.py",{'dir':'ForceField'}),
  ("python","test_list.py",{'dir':'DistGeom'}),
  ("python","test_list.py",{'dir':'ChemicalFeatures'}),
  ("python","test_list.py",{'dir':'ChemicalFeatures/Wrap'}),
  #("python","test_list.py",{'dir':'PgSQL/RDLib'}),
  ("python","test_list.py",{'dir':'SimDivPickers/Wrap'}),
  ("python","test_list.py",{'dir':'ML'}),
  ]

longTests = [
  ("python","test_list.py",{'dir':'GraphMol'}),
  ]

if __name__=='__main__':
  import sys
  import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
