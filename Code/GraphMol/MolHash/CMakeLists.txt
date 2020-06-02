
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLHASH_BUILD)
rdkit_library(MolHash
              hashfunctions.cpp normalize.cpp
              LINK_LIBRARIES Depictor Descriptors GraphMol)

rdkit_headers(MolHash.h nmmolhash.h           
              DEST GraphMol/MolHash)

rdkit_catch_test(molHashCatchTest catch_tests.cpp LINK_LIBRARIES MolHash)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

