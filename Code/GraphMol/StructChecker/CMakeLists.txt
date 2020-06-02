
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_STRUCTCHECKER_BUILD)
rdkit_library(StructChecker
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp Utilites.cpp
              Pattern.cpp Stereo.cpp AtomSymbolMatch.cpp StripSmallFragments.cpp ReCharge.cpp Tautomer.cpp
              LINK_LIBRARIES Depictor SubstructMatch FileParsers ChemTransforms Descriptors
                             GraphMol RDGeneral )

rdkit_headers(StructChecker.h Utilites.h Pattern.h Stereo.h StripSmallFragments.h ReCharge.h Tautomer.h
              DEST GraphMol/StructChecker)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker SmilesParse )
