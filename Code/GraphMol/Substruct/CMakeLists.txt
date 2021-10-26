
rdkit_library(SubstructMatch 
              SubstructMatch.cpp SubstructUtils.cpp Generics.cpp
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(SubstructMatch PRIVATE RDKIT_SUBSTRUCTMATCH_BUILD)

rdkit_headers(SubstructMatch.h
              SubstructUtils.h Generics.h DEST GraphMol/Substruct)

rdkit_test(testSubstructMatch test1.cpp LINK_LIBRARIES FileParsers SmilesParse SubstructMatch)

rdkit_catch_test(substructTestCatch ../catch_main.cpp catch_tests.cpp LINK_LIBRARIES FileParsers SmilesParse SubstructMatch)
rdkit_catch_test(substructTestGenerics ../catch_main.cpp generic_tests.cpp LINK_LIBRARIES FileParsers SmilesParse SubstructMatch)
