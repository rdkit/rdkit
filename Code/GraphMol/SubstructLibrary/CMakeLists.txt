if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making SubstructLibrary without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

if(RDK_BUILD_THREADSAFE_SSS)

remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
	      PatternFactory.cpp
              LINK_LIBRARIES Fingerprints SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral DataStructs ${RDKit_SERIALIZATION_LIBS})

rdkit_headers(SubstructLibrary.h
              SubstructLibrarySerialization.h
              PatternFactory.h
              DEST GraphMol/SubstructLibrary)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs
           ${RDKit_SERIALIZATION_LIBS})

endif(RDK_BUILD_THREADSAFE_SSS)
