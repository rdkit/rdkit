
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SUBSTRUCTMATCH_BUILD)
rdkit_library(SubstructMatch 
              SubstructMatch.cpp SubstructUtils.cpp
              LINK_LIBRARIES GraphMol
                ${RDKit_THREAD_LIBS} )

rdkit_headers(SubstructMatch.h
              SubstructUtils.h DEST GraphMol/Substruct)

rdkit_test(testSubstructMatch test1.cpp LINK_LIBRARIES  FileParsers SmilesParse SubstructMatch
GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS} )

rdkit_catch_test(substructTestCatch catch_tests.cpp LINK_LIBRARIES FileParsers SmilesParse 
SubstructMatch GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
