add_library(FileParsers MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp 
  ForwardSDMolSupplier.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp 
  SmilesWriter.cpp SDWriter.cpp TDTWriter.cpp TDTMolSupplier.cpp 
  TplFileParser.cpp TplFileWriter.cpp Mol2FileParser.cpp )

add_executable(fileParsersTest1 test1.cpp)
target_link_libraries(fileParsersTest1  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(fileParsersTest1 ${EXECUTABLE_OUTPUT_PATH}/fileParsersTest1 )

add_executable(testMolSupplier testMolSupplier.cpp)
target_link_libraries(testMolSupplier  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(testMolSupplier ${EXECUTABLE_OUTPUT_PATH}/testMolSupplier )

add_executable(testMolWriter testMolWriter.cpp)
target_link_libraries(testMolWriter  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(testMolWriter ${EXECUTABLE_OUTPUT_PATH}/testMolWriter )

add_executable(testTplParser testTpls.cpp)
target_link_libraries(testTplParser  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(testTplParser ${EXECUTABLE_OUTPUT_PATH}/testTplParser )

add_executable(testMol2ToMol testMol2ToMol.cpp)
target_link_libraries(testMol2ToMol  FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(testMol2ToMol ${EXECUTABLE_OUTPUT_PATH}/testMol2ToMol )

