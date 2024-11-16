
rdkit_library(SynthonSpaceSearch
        SynthonSpaceSearch_details.cpp SynthonSpace.cpp SynthonSet.cpp Synthon.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms Fingerprints SubstructMatch GraphMol)
target_compile_definitions(SynthonSpaceSearch PRIVATE RDKIT_SYNTHONSPACESEARCH_BUILD)

rdkit_headers(SynthonSpace.h SynthonSet.h Synthon.h SubstructureResults.h
        DEST GraphMol/SynthonSpaceSearch)

rdkit_catch_test(testSynthonSpaceSearch synthon_space_search_catch.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
