rdkit_library(FragCatalog
              FragCatalogUtils.cpp FragCatGenerator.cpp FragCatalogEntry.cpp 
              FragCatParams.cpp FragFPGenerator.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse Catalogs GraphMol RDGeometryLib RDGeneral )

rdkit_test(testFragCatalog test1.cpp LINK_LIBRARIES FragCatalog Subgraphs SubstructMatch SmilesParse Catalogs FileParsers GraphMol DataStructs RDGeometryLib RDGeneral  )

add_subdirectory(Wrap)

