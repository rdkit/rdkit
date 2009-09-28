rdkit_library(FragCatalog
              FragCatalogUtils.cpp FragCatGenerator.cpp FragCatalogEntry.cpp 
              FragCatParams.cpp FragFPGenerator.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse Catalogs)

rdkit_test(testFragCatalog test1.cpp LINK_LIBRARIES FileParsers FragCatalog)

add_subdirectory(Wrap)

