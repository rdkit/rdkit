if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(SUBSTRUCTLIBRARY_USE_BOOST_SERIALIZATION "yes")
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making SubstructLibrary without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

configure_file(SubstructLibraryDefs.h.in SubstructLibraryDefs.h @ONLY)

if(RDK_BUILD_THREADSAFE_SSS)

remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
              LINK_LIBRARIES Fingerprints SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral DataStructs ${RDKit_SERIALIZATION_LIBS})

rdkit_headers(SubstructLibrary.h
              SubstructLibraryDefs.h
              DEST GraphMol/SubstructLibrary)

add_subdirectory(Wrap)

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs
           ${RDKit_SERIALIZATION_LIBRARY})
endif(RDK_BUILD_THREADSAFE_SSS)
