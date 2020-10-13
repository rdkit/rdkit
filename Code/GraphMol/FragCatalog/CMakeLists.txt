
rdkit_library(FragCatalog
              FragCatalogUtils.cpp FragCatGenerator.cpp FragCatalogEntry.cpp 
              FragCatParams.cpp FragFPGenerator.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse Catalogs GraphMol RDGeneral )
target_compile_definitions(FragCatalog PRIVATE RDKIT_FRAGCATALOG_BUILD)

rdkit_headers(FragCatalogEntry.h
              FragCatalogUtils.h
              FragCatGenerator.h
              FragCatParams.h
              FragFPGenerator.h DEST GraphMol/FragCatalog)

rdkit_test(testFragCatalog test1.cpp LINK_LIBRARIES FragCatalog FileParsers )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
