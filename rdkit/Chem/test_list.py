
tests=[
  ("python","UnitTestChem.py",{}),
  ("python","UnitTestChemv2.py",{}),
  ("python","UnitTestChemAtom.py",{}),
  ("python","UnitTestChemBond.py",{}),
  ("python","UnitTestChemSmarts.py",{}),
  ("python","UnitTestFragmentDescriptors.py",{}),
  ("python","UnitTestGraphDescriptors.2.py",{}),
  ("python","UnitTestLipinski.py",{}),
  ("python","UnitTestOldBugs.py",{}),
  ("python","UnitTestSATIS.py",{}),
  ("python","UnitTestSmiles.py",{}),
  ("python","UnitTestSuppliers.py",{}),
  ("python","UnitTestSurf.py",{}),
  ("python","FragmentMatcher.py",{}),
  ("python","MACCSkeys.py",{}),
  ("python","Descriptors.py",{}),
  ("python","UnitTestCatalog.py",{}),
  ("python","TemplateAlign.py",{}),
  ("python","Recap.py",{}),
  ("python","UnitTestDescriptors.py",{}),
  ("python","AllChem.py",{}),
  ("python","test_list.py",{'dir':'AtomPairs'}),
  ("python","test_list.py",{'dir':'ChemUtils'}),
  ("python","test_list.py",{'dir':'EState'}),
  ("python","test_list.py",{'dir':'FeatMaps'}),
  ("python","test_list.py",{'dir':'Fingerprints'}),
  ("python","test_list.py",{'dir':'Pharm2D'}),
  ("python","test_list.py",{'dir':'Pharm3D'}),
  ("python","test_list.py",{'dir':'Subshape'}),
  ("python","test_list.py",{'dir':'Suppliers'}),
  ]



longTests=[
  ("python","UnitTestArom.py",{}),
  ("python","UnitTestCrippen.py",{}),
  ("python","UnitTestGraphDescriptors.2.py -l",{}),
  ("python","UnitTestSurf.py -l",{}),
  ]

if __name__=='__main__':
  import sys
  import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
