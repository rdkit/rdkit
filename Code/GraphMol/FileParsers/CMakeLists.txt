add_library(FileParsers MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
  ForwardSDMolSupplier.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp 
  SmilesWriter.cpp SDWriter.cpp TDTWriter.cpp TDTMolSupplier.cpp 
  TplFileParser.cpp TplFileWriter.cpp Mol2FileParser.cpp )

rdkit_test(fileParsersTest1 test1.cpp LINK_LIBRARIES 
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
rdkit_test(testMolSupplier testMolSupplier.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
rdkit_test(testTplParser testTpls.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
rdkit_test(testMol2ToMol testMol2ToMol.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
