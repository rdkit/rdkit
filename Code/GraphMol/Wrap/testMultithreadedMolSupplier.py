import os, sys
import unittest, doctest
from rdkit import RDConfig, rdBase
from rdkit import Chem


class TestCase(unittest.TestCase): 
  def testMultiSmiMolSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data', 'first_200.tpsa.csv')
    #fileN = "../FileParsers/test_data/first_200.tpsa.csv"
    smiSup = Chem.MultithreadedSmilesMolSupplier(fileN, ",", 0, - 1) 
    i = 0
    while not smiSup.atEnd():
      mol = next(smiSup)
      if(mol):
        i += 1
    self.assertTrue(i == 200)
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data', 'fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"
    smiSup = Chem.MultithreadedSmilesMolSupplier(fileN, delimiter = ",", smilesColumn = 1, nameColumn = 0, titleLine = 0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    confusedNames = []
    i = 0
    while not smiSup.atEnd():
      mol = next(smiSup)
      confusedNames.append(mol.GetProp("_Name"))
      i += 1
    self.assertTrue(i == 10)
    self.assertTrue(sorted(confusedNames) == sorted(names))


  def testMultiSDMolSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data', 'NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    sdSup = Chem.MultithreadedSDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213", "220", "229", "256"]
    confusedMolNames = []
    i = 0
    for mol in sdSup:
      if(mol):
        confusedMolNames.append(mol.GetProp("_Name"))
        i += 1 
    self.assertTrue(len(molNames) == i)
    self.assertTrue(sorted(confusedMolNames) == sorted(molNames))


if __name__ == '__main__':
  print("Testing MultithreadedMolSupplier (Smiles and SD)")
  unittest.main()
