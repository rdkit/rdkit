add_library(FileParsers MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
  ForwardSDMolSupplier.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp 
  SmilesWriter.cpp SDWriter.cpp )

rdkit_test(fileParsersTest1 test1.cpp LINK_LIBRARIES 
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )
rdkit_test(testMolSupplier testMolSupplier.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )
rdkit_test(testMolWriter testMolWriter.cpp LINK_LIBRARIES
  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )
