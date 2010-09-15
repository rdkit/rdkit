rdkit_library(Catalogs 
              Catalog.cpp CatalogParams.cpp CatalogEntry.cpp
              LINK_LIBRARIES RDGeneral)

rdkit_headers(CatalogEntry.h
              Catalog.h
              CatalogParams.h DEST Catalogs)


