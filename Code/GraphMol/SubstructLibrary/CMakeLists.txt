if(RDK_BUILD_THREADSAFE_SSS)
rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
              LINK_LIBRARIES Fingerprints SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral DataStructs )

rdkit_headers(SubstructLibrary.h
              DEST GraphMol/SubstructLibrary)

add_subdirectory(Wrap)

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs
           ${Boost_SERIALIZATION_LIBRARY})
endif(RDK_BUILD_THREADSAFE_SSS)
