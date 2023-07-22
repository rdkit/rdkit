
rdkit_library(RascalMCES
        rascal_mces.cpp lap_a_la_scipy.cpp
        PartitionSet.cpp RascalResult.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms SubstructMatch GraphMol)
target_compile_definitions(RascalMCES PRIVATE RDKIT_RASCALMCES_BUILD)

rdkit_headers(rascal_mces.h PartitionSet.h RascalOptions.h RascalResult.h
        DEST GraphMol/RascalMCES)

rdkit_test(testRascalMCES mces_catch.cpp LINK_LIBRARIES RascalMCES)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
#add_subdirectory(Test)
