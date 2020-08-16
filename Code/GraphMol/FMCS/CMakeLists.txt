
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_FMCS_BUILD)
rdkit_library(FMCS
              FMCS.cpp Seed.cpp MaximumCommonSubgraph.cpp SubstructMatchCustom.cpp
              LINK_LIBRARIES Depictor SmilesParse FileParsers ChemTransforms SubstructMatch GraphMol )

rdkit_headers(FMCS.h
              Graph.h
              DEST GraphMol/FMCS)

rdkit_test(testFMCS testFMCS_Unit.cpp LINK_LIBRARIES
FMCS )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
add_subdirectory(Test)
