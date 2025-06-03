
rdkit_library(EnumerateStereoisomers
        EnumerateStereoisomers.cpp Flippers.cpp
        LINK_LIBRARIES SmilesParse FileParsers GraphMol CIPLabeler DistGeomHelpers DistGeometry)
target_compile_definitions(EnumerateStereoisomers PRIVATE RDKIT_ENUMERATESTEREOISOMERS_BUILD)

rdkit_headers(EnumerateStereoisomers.h Flippers.h DEST GraphMol/EnumerateStereoisomers)

rdkit_catch_test(testEnumerateStereoisomers catch_tests.cpp
        LINK_LIBRARIES EnumerateStereoisomers)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
