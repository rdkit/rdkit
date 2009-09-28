rdkit_library(FileParsers
              Mol2FileParser.cpp  
              MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
              ForwardSDMolSupplier.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp 
              SmilesWriter.cpp SDWriter.cpp TDTMolSupplier.cpp TDTWriter.cpp
              TplFileParser.cpp TplFileWriter.cpp
              LINK_LIBRARIES SmilesParse)

rdkit_test(fileParsersTest1 test1.cpp 
           LINK_LIBRARIES FileParsers Depictor SubstructMatch)

rdkit_test(testMolSupplier testMolSupplier.cpp 
           LINK_LIBRARIES FileParsers Depictor SubstructMatch)

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES FileParsers)

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers)

