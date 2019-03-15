



if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
    add_definitions(-DRDK_USE_BOOST_SERIALIZATION)
else()
    message("== Making SubstructLibrary without boost Serialization support")
    remove_definitions(-DRDK_USE_BOOST_SERIALIZATION)
    set(RDKit_SERIALIZATION_LIBS )
endif()

if(RDK_BUILD_THREADSAFE_SSS)

find_package(Boost REQUIRED COMPONENTS system thread)

remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD -DBOOST_SYSTEM_NO_DEPRECATED)
rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
              EditableSubstructLibrary
              EditableSubstructLibraryTrustedSmilesWithPattern.cpp
              ChunkedHitlist
              LINK_LIBRARIES Fingerprints SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral DataStructs ${RDKit_SERIALIZATION_LIBS})

rdkit_headers(SubstructLibrary.h
              SubstructLibrarySerialization.h
              DEST GraphMol/SubstructLibrary)

add_subdirectory(Wrap)

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs
           ${Boost_SERIALIZATION_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})

rdkit_test(editableSubstructLibraryTest editableSubstructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs
           ${Boost_SERIALIZATION_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})

endif(RDK_BUILD_THREADSAFE_SSS)
