
tests=[
  ("testExecs/test1.exe","",{}),
  ("testExecs/testMolSupplier.exe","",{}),
  ("testExecs/testMolWriter.exe","",{}),
  ("testExecs/testTplParser.exe","",{}),
  ("testExecs/testMol2ToMol.exe","",{}),
  ]



longTests=[]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
