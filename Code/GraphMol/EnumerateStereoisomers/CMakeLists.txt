
rdkit_library(EnumerateStereoisomers
        EnumerateStereoisomers.cpp Flippers.cpp RingSystemFilter.cpp
        LINK_LIBRARIES SmilesParse FileParsers GraphMol CIPLabeler DistGeomHelpers DistGeometry)
target_compile_definitions(EnumerateStereoisomers PRIVATE RDKIT_ENUMERATESTEREOISOMERS_BUILD)

rdkit_headers(EnumerateStereoisomers.h Flippers.h DEST GraphMol/EnumerateStereoisomers)

rdkit_catch_test(testEnumerateStereoisomers catchEnumerateStereoisomers.cpp RingSystemFilter.cpp
        LINK_LIBRARIES EnumerateStereoisomers)

if (RDK_BUILD_BOOST_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()
