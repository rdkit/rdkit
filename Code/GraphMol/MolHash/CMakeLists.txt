
rdkit_library(MolHash
              hashfunctions.cpp normalize.cpp
              LINK_LIBRARIES SubstructMatch Depictor GraphMol)
target_compile_definitions(MolHash PRIVATE RDKIT_MOLHASH_BUILD)

rdkit_headers(MolHash.h nmmolhash.h           
              DEST GraphMol/MolHash)

rdkit_catch_test(molHashCatchTest catch_tests.cpp  LINK_LIBRARIES SmilesParse FileParsers MolHash)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

