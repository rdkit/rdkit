
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_STRUCTCHECKER_BUILD)
rdkit_library(StructChecker
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp Utilites.cpp
              Pattern.cpp Stereo.cpp AtomSymbolMatch.cpp StripSmallFragments.cpp ReCharge.cpp Tautomer.cpp
              LINK_LIBRARIES Depictor SubstructMatch FileParsers ChemTransforms Descriptors
                             GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS})

rdkit_headers(StructChecker.h Utilites.h Pattern.h Stereo.h StripSmallFragments.h ReCharge.h Tautomer.h
              DEST GraphMol/StructChecker)

add_subdirectory(Wrap)

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker ChemTransforms Depictor Descriptors
                          SubstructMatch FileParsers SmilesParse 
                          GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
