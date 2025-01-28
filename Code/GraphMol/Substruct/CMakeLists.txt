
rdkit_library(SubstructMatch 
              SubstructMatch.cpp SubstructUtils.cpp 
              LINK_LIBRARIES GenericGroups GraphMol RDGeneral)
target_compile_definitions(SubstructMatch PRIVATE RDKIT_SUBSTRUCTMATCH_BUILD)

rdkit_headers(SubstructMatch.h
              SubstructUtils.h DEST GraphMol/Substruct)

rdkit_test(testSubstructMatch test1.cpp LINK_LIBRARIES FileParsers SmilesParse SubstructMatch)

rdkit_catch_test(substructTestCatch catch_tests.cpp LINK_LIBRARIES FileParsers SmilesParse SubstructMatch)
