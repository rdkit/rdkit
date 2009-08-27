
tests=[
  ("python","UnitTestChem.py",{}),
  ("python","UnitTestChemv2.py",{}),
  ("python","UnitTestChemAtom.py",{}),
  ("python","UnitTestChemBond.py",{}),
  ("python","UnitTestChemSmarts.py",{}),
  ("python","UnitTestOldBugs.py",{}),
  ("python","UnitTestSmiles.py",{}),
  ("python","UnitTestSuppliers.py",{}),
  ("python","AllChem.py",{}),
  ("python","PropertyMol.py",{}),
  ]



longTests=[
  ("python","UnitTestArom.py",{}),
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
