if(RDK_BUILD_THREADSAFE_SSS)
rdkit_library(SynthonSpaceSearch
        SynthonSpaceSearch_details.cpp SynthonSpace.cpp SynthonSet.cpp Synthon.cpp
        SynthonSpaceSearcher.cpp SynthonSpaceSubstructureSearcher.cpp SynthonSpaceFingerprintSearcher.cpp
        SynthonSpaceRascalSearcher.cpp SynthonSpaceShapeSearcher.cpp
        SynthonSpaceHitSet.cpp SearchResults.cpp SynthonShapeInput.cpp MemoryMappedFileReader.cpp
        LINK_LIBRARIES SmilesParse FileParsers ChemTransforms Fingerprints SubstructMatch GraphMol RascalMCES
        GeneralizedSubstruct DistGeomHelpers DistGeometry SimDivPickers Descriptors
        EnumerateStereoisomers GaussianShape)
target_compile_definitions(SynthonSpaceSearch PRIVATE RDKIT_SYNTHONSPACESEARCH_BUILD)

rdkit_headers(SynthonSpace.h SynthonSet.h Synthon.h SearchResults.h SynthonShapeInput.h
        SynthonSpaceSearcher.h SynthonSpaceSubstructureSearcher.h SynthonSpaceFingerprintSearcher.h
        SynthonSpaceRascalSearcher.h SynthonSpaceShapeSearcher.h ProgressBar.h DEST GraphMol/SynthonSpaceSearch)

rdkit_catch_test(testSynthonSpaceSubstructureSearch substructure_search_catch_tests.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)
rdkit_catch_test(testSynthonSpaceFingerprintSearch fingerprint_search_catch_tests.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)
rdkit_catch_test(testSynthonSpaceRascalSearch rascal_search_catch_tests.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)
rdkit_catch_test(testSynthonSpaceShapeSearch shape_search_catch_tests.cpp
        LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)

if (RDK_BUILD_LONG_RUNNING_TESTS)
    rdkit_catch_test(testSynthonSpaceLONG_TEST LONG_TEST_catch_tests.cpp
            LINK_LIBRARIES SubstructMatch SubstructLibrary SynthonSpaceSearch)
endif ()

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
else()
message("The SynthonSpaceSearch module requires threads, so not building.")
endif()