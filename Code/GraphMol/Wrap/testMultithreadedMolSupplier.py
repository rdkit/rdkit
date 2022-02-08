import os, sys, unittest, doctest
import gzip
from rdkit import RDConfig, rdBase
from rdkit import Chem
from rdkit import __version__
import sys

class TestCase(unittest.TestCase):
    def testMultiSmiMolSupplier(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol',
                             'FileParsers', 'test_data', 'first_200.tpsa.csv')
        # fileN = "../FileParsers/test_data/first_200.tpsa.csv"
        smiSup = Chem.MultithreadedSmilesMolSupplier(fileN, ",", 0, - 1)
        i = 0
        while not smiSup.atEnd():
            mol = next(smiSup)
            if(mol):
                i += 1
        self.assertTrue(i == 200)
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol',
                             'FileParsers', 'test_data', 'fewSmi.csv')
        # fileN = "../FileParsers/test_data/fewSmi.csv"
        smiSup = Chem.MultithreadedSmilesMolSupplier(
            fileN, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
        names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
        props = ["34.14", "25.78", "106.51", "82.78", "60.16",
            "87.74", "37.38", "77.28", "65.18", "0.00"]
        confusedNames = []
        confusedProps = []
        i = 0
        for mol in smiSup:
            if mol is not None:
                self.assertTrue(mol.HasProp("_Name"))
                self.assertTrue(mol.HasProp("Column_2"))
                prop = mol.GetProp("Column_2")
                name = mol.GetProp("_Name")
                confusedProps.append(prop)
                confusedNames.append(name)
                i += 1
        self.assertTrue(i == 10)
        self.assertTrue(sorted(confusedNames) == sorted(names))
        self.assertTrue(sorted(confusedProps) == sorted(props))

        # context manager
        confusedNames = []
        confusedProps = []
        i = 0
        with Chem.MultithreadedSmilesMolSupplier(fileN,delimiter=",", smilesColumn=1, 
                  nameColumn=0, titleLine=0) as smiSup:
            for mol in smiSup:
                if mol is not None:
                    self.assertTrue(mol.HasProp("_Name"))
                    self.assertTrue(mol.HasProp("Column_2"))
                    prop = mol.GetProp("Column_2")
                    name = mol.GetProp("_Name")
                    confusedProps.append(prop)
                    confusedNames.append(name)
                    i += 1
            self.assertTrue(i == 10)
            self.assertTrue(sorted(confusedNames) == sorted(names))
            self.assertTrue(sorted(confusedProps) == sorted(props))


    def testMultiSDMolSupplier(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol',
                             'FileParsers', 'test_data', 'NCI_aids_few.sdf')
        # fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
        sdSup = Chem.MultithreadedSDMolSupplier(fileN)
        molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]
        confusedMolNames = []
        i = 0
        for mol in sdSup:
            if mol is not None:
                confusedMolNames.append(mol.GetProp("_Name"))
                i += 1
        self.assertTrue(len(molNames) == i)
        self.assertTrue(sorted(confusedMolNames) == sorted(molNames))

        # context manager
        confusedMolNames = []
        i = 0
        with Chem.MultithreadedSDMolSupplier(fileN) as sdSup:
            for mol in sdSup:
                if mol is not None:
                    confusedMolNames.append(mol.GetProp("_Name"))
                    i += 1
        self.assertTrue(len(molNames) == i)
        self.assertTrue(sorted(confusedMolNames) == sorted(molNames))





    # NOTE these are disabled until we rewrite the code to construct a 
    #      MultithreadedSDMolSupplier from a python stream
    @unittest.skip("Skipping construction from stream")
    def testMultiSDMolSupplierFromStream(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol',
                             'FileParsers', 'test_data', 'NCI_aids_few.sdf')
        molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]
        # try opening with streambuf
        inf = open(fileN,'rb')
        if(inf):
          gSup = Chem.SDMolSupplierFromStream(inf)
          confusedMolNames = []
          i = 0
          for mol in gSup:
            # print("!!",i,file=sys.stderr);sys.stderr.flush()
            if(mol):
              confusedMolNames.append(mol.GetProp("_Name"))
              i += 1
          self.assertTrue(len(molNames) == i)
          self.assertTrue(sorted(confusedMolNames) == sorted(molNames))
        #   print("done!",file=sys.stderr);sys.stderr.flush()
        # try opening with streambuf
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz') 
        # try opening with gzip
        inf = gzip.open(fileN)
        if(inf):
          gSup = Chem.SDMolSupplierFromStream(inf)
          confusedMolNames = []
          i = 0
          for mol in gSup:
            # print("!",i,file=sys.stderr);sys.stderr.flush()
            if(mol):
              confusedMolNames.append(mol.GetProp("_Name"))
              i += 1
          self.assertTrue(len(molNames) == i)
          self.assertTrue(sorted(confusedMolNames) == sorted(molNames))
           
          

  
if __name__ == '__main__':
    print("Testing Smiles and SD MultithreadedMolSupplier")
    unittest.main()
