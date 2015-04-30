
tests=[
  ("python","UnitTestChem.py",{}),
  ("python","UnitTestChemv2.py",{}),
  ("python","UnitTestChemAtom.py",{}),
  ("python","UnitTestChemBond.py",{}),
  ("python","UnitTestChemSmarts.py",{}),
  ("python","UnitTestFragmentDescriptors.py",{}),
  ("python","UnitTestGraphDescriptors.2.py",{}),
  ("python","UnitTestLipinski.py",{}),
  ("python","MCS.py",{}),
  ("python","UnitTestMCS.py",{}),
  ("python","UnitTestOldBugs.py",{}),
  ("python","UnitTestSATIS.py",{}),
  ("python","UnitTestSmiles.py",{}),
  ("python","UnitTestSuppliers.py",{}),
  ("python","UnitTestSurf.py",{}),
  ("python","UnitTestMol3D.py",{}),
  ("python","FragmentMatcher.py",{}),
  ("python","MACCSkeys.py",{}),
  ("python","Descriptors.py",{}),
  ("python","UnitTestCatalog.py",{}),
  ("python","TemplateAlign.py",{}),
  ("python","Recap.py",{}),
  ("python","BRICS.py",{}),
  ("python","UnitTestDescriptors.py",{}),
  ("python","AllChem.py",{}),
  ("python","PropertyMol.py",{}),
  ("python","UnitTestInchi.py",{}),
  ("python","SaltRemover.py",{}),
  ("python","UnitTestFunctionalGroups.py",{}),
  ("python","UnitTestCrippen.py",{}),
  ("python","__init__.py",{}),
  ("python","PandasTools.py",{}),
  ("python","test_list.py",{'dir':'AtomPairs'}),
  ("python","test_list.py",{'dir':'ChemUtils'}),
  ("python","test_list.py",{'dir':'EState'}),
  ("python","test_list.py",{'dir':'FeatMaps'}),
  ("python","test_list.py",{'dir':'Fingerprints'}),
  ("python","test_list.py",{'dir':'Pharm2D'}),
  ("python","test_list.py",{'dir':'Pharm3D'}),
  #("python","test_list.py",{'dir':'Subshape'}),
  ("python","test_list.py",{'dir':'Suppliers'}),
  ("python","test_list.py",{'dir':'Scaffolds'}),
  ("python","test_list.py",{'dir':'Draw'}),
  ("python","test_list.py",{'dir':'Fraggle'}),
  ("python","test_list.py",{'dir':'SimpleEnum'}),
  ]



longTests=[
  ("python","UnitTestArom.py",{}),
  ("python","UnitTestGraphDescriptors.2.py -l",{}),
  ("python","UnitTestSurf.py -l",{}),
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
