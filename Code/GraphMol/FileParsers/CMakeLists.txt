rdkit_library(FileParsers
              Mol2FileParser.cpp  
              MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
              ForwardSDMolSupplier.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp 
              SmilesWriter.cpp SDWriter.cpp TDTMolSupplier.cpp TDTWriter.cpp
              TplFileParser.cpp TplFileWriter.cpp)
              

rdkit_test(fileParsersTest1 test1.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMolstatic RDGeneral RDGeometryLib )

rdkit_test(testMolSupplier testMolSupplier.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMolstatic RDGeneral RDGeometryLib )

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMolstatic RDGeneral RDGeometryLib )

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMolstatic RDGeneral RDGeometryLib )

rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMolstatic RDGeneral RDGeometryLib )
