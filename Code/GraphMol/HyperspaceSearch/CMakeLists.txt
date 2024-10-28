
rdkit_library(HyperspaceSearch
        HyperspaceSubstructure_details.cpp Hyperspace.cpp ReactionSet.cpp Reagent.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms Fingerprints SubstructMatch GraphMol)
target_compile_definitions(HyperspaceSearch PRIVATE RDKIT_HYPERSPACE_SEARCH_BUILD)

rdkit_headers(HyperspaceSubstructure_details.h Hyperspace.h ReactionSet.h Reagent.h SubstructureResults.h
        DEST GraphMol/HyperspaceSearch)

rdkit_catch_test(testHyperspaceSearch hyperspace_search_catch.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary HyperspaceSearch)

add_executable(convertSynthonSpace convert_synthon_space.cpp)
target_link_libraries(convertSynthonSpace SubstructMatch SubstructLibrary HyperspaceSearch)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
