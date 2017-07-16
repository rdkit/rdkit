rdkit_library(SubstructLibrary
              SubstructLibrary.cpp
              LINK_LIBRARIES Fingerprints SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral)

rdkit_headers(SubstructLibrary.h
              DEST GraphMol/SubstructLibrary)

add_subdirectory(Wrap)

rdkit_test(substructLibraryTest substructLibraryTest.cpp
           LINK_LIBRARIES SubstructLibrary Fingerprints
           SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral
           ${Boost_SERIALIZATION_LIBRARY})
