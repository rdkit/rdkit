
rdkit_library(RascalMCES
        RascalMCES.cpp RascalCluster.cpp RascalButinaCluster.cpp
        lap_a_la_scipy.cpp PartitionSet.cpp RascalResult.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms SubstructMatch GraphMol)
target_compile_definitions(RascalMCES PRIVATE RDKIT_RASCALMCES_BUILD)

rdkit_headers(RascalMCES.h RascalOptions.h RascalClusterOptions.h RascalResult.h
        DEST GraphMol/RascalMCES)

rdkit_catch_test(testRascalMCES mces_catch.cpp LINK_LIBRARIES RascalMCES)
rdkit_catch_test(testRascalCluster mces_cluster_catch.cpp LINK_LIBRARIES RascalMCES)
if (RDK_BUILD_LONG_RUNNING_TESTS)
    rdkit_catch_test(testRascalMCESLONG_TEST LONG_TEST_mces_catch.cpp LINK_LIBRARIES RascalMCES)
endif ()

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
