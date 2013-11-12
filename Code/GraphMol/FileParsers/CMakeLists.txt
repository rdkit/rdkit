rdkit_library(FileParsers
              Mol2FileParser.cpp  
              MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
              ForwardSDMolSupplier.cpp SDMolSupplier.cpp SDWriter.cpp
              SmilesMolSupplier.cpp SmilesWriter.cpp
              TDTMolSupplier.cpp TDTWriter.cpp
              TplFileParser.cpp TplFileWriter.cpp
              PDBParser.cpp PDBWriter.cpp PDBSupplier.cpp ProximityBonds.cpp
              LINK_LIBRARIES SmilesParse GraphMol)
              
rdkit_headers(FileParsers.h
              FileParserUtils.h
              MolFileStereochem.h
              MolSupplier.h
              MolWriters.h DEST GraphMol/FileParsers)

rdkit_test(fileParsersTest1 test1.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMol RDGeneral RDGeometryLib )

rdkit_test(testMolSupplier testMolSupplier.cpp 
           LINK_LIBRARIES FileParsers SmilesParse Depictor SubstructMatch GraphMol RDGeneral RDGeometryLib )

rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )
