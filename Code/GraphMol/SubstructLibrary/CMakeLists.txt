if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making SubstructLibrary without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
	      PatternFactory.cpp
              LINK_LIBRARIES  TautomerQuery MolStandardize Fingerprints SubstructMatch SmilesParse
              GraphMol Catalogs DataStructs RDGeneral ${RDKit_SERIALIZATION_LIBS})
target_compile_definitions(SubstructLibrary PRIVATE RDKIT_SUBSTRUCTLIBRARY_BUILD)

rdkit_headers(SubstructLibrary.h
              SubstructLibrarySerialization.h
              PatternFactory.h
              DEST GraphMol/SubstructLibrary)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES TautomerQuery MolStandardize SubstructLibrary FileParsers )

rdkit_catch_test(substructLibraryCatchTest catch_tests.cpp LINK_LIBRARIES SubstructLibrary FileParsers )
