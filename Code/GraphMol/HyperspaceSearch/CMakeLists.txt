
rdkit_library(HyperspaceSearch
        HyperspaceSubstructureSearch.cpp Hyperspace.cpp Reagent.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms Fingerprints SubstructMatch GraphMol)
target_compile_definitions(HyperspaceSearch PRIVATE RDKIT_HYPERSPACE_SEARCH_BUILD)

rdkit_headers(HyperspaceSubstructureSearch.h Hyperspace.h Reagent.h
        DEST GraphMol/HyperspaceSearch)

rdkit_catch_test(testHyperspaceSearch hyperspace_search_catch.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary HyperspaceSearch)

#if (RDK_BUILD_PYTHON_WRAPPERS)
#    add_subdirectory(Wrap)
#endif ()
