import sys

tests=[
  ("testExecs/itertest.exe","",{}),
  ("testExecs/MolOpsTest.exe","",{}),
  ("testExecs/testCanon.exe","C1OCCC1 C1CCOC1",{}),
  ("testExecs/testPickler.exe","",{}),
  ("testExecs/test1.exe","",{}),
  ("testExecs/testChirality.exe","",{}),

  ("python","test_list.py",{'dir':'Depictor'}),
  ("python","test_list.py",{'dir':'FileParsers'}),
  ("python","test_list.py",{'dir':'SmilesParse'}),
  ("python","test_list.py",{'dir':'Substruct'}),
  ("python","test_list.py",{'dir':'Subgraphs'}),
  ("python","test_list.py",{'dir':'FragCatalog'}),
  ("python","test_list.py",{'dir':'Fingerprints'}),
  ("python","test_list.py",{'dir':'MolTransforms'}),

  ("python","test_list.py",{'dir':'Wrap'}),
  ("python","test_list.py",{'dir':'Depictor/Wrap'}),
  ("python","test_list.py",{'dir':'FragCatalog/Wrap'}),
  ("python","test_list.py",{'dir':'PartialCharges/Wrap'}),

  ("python","test_list.py",{'dir':'ForceFieldHelpers'}),
  ("python","test_list.py",{'dir':'DistGeomHelpers'}),
  ("python","test_list.py",{'dir':'Descriptors'}),
  ("python","test_list.py",{'dir':'Descriptors/Wrap'}),
  ("python","test_list.py",{'dir':'MolChemicalFeatures'}),
  ("python","test_list.py",{'dir':'MolAlign'}),
  ("python","test_list.py",{'dir':'ShapeHelpers'}),
  ("python","test_list.py",{'dir':'ChemTransforms'}),

  ("python","test_list.py",{'dir':'MolCatalog'}),
  ("python","test_list.py",{'dir':'MolCatalog/Wrap'}),

  ("python","test_list.py",{'dir':'ChemReactions'}),

  ("python","test_list.py",{'dir':'SLNParse'}),
  ("python","test_list.py",{'dir':'SLNParse/Wrap'}),
  
  ]

if sys.platform != 'win32':
  tests.extend([
    ("testExecs/cptest.exe","",{}),
    ("testExecs/querytest.exe","",{}),
    ])



longTests=[
  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
