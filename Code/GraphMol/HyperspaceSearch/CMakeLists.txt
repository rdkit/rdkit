
rdkit_library(HyperspaceSearch
        HyperspaceSubstructureSearch.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms SubstructMatch GraphMol)
target_compile_definitions(HyperspaceSearch PRIVATE RDKIT_HYPERSPACE_SEARCH_BUILD)

rdkit_headers(HyperspaceSubstructureSearch.h
        DEST GraphMol/HyperspaceSearch)

rdkit_catch_test(testHyperspaceSearch hyperspace_search_catch.cpp LINK_LIBRARIES HyperspaceSearch)

#if (RDK_BUILD_PYTHON_WRAPPERS)
#    add_subdirectory(Wrap)
#endif ()
