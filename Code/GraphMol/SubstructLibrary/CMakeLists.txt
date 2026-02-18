if(NOT RDK_USE_BOOST_SERIALIZATION)
    message("== Making SubstructLibrary without boost Serialization support")
endif()

rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
	      PatternFactory.cpp
              LINK_LIBRARIES  GeneralizedSubstruct TautomerQuery MolStandardize Fingerprints SubstructMatch SmilesParse
              GraphMol Catalogs DataStructs RDGeneral)
target_compile_definitions(SubstructLibrary PRIVATE RDKIT_SUBSTRUCTLIBRARY_BUILD)

rdkit_headers(SubstructLibrary.h
              SubstructLibrarySerialization.h
              PatternFactory.h
              DEST GraphMol/SubstructLibrary)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_catch_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES TautomerQuery MolStandardize SubstructLibrary FileParsers )

rdkit_catch_test(substructLibraryCatchTest catch_tests.cpp LINK_LIBRARIES SubstructLibrary FileParsers )
