
rdkit_library(Catalogs 
              Catalog.cpp CatalogParams.cpp CatalogEntry.cpp
              LINK_LIBRARIES RDGeneral)
target_compile_definitions(Catalogs PRIVATE RDKIT_CATALOGS_BUILD)

rdkit_headers(CatalogEntry.h
              Catalog.h
              CatalogParams.h DEST Catalogs)


