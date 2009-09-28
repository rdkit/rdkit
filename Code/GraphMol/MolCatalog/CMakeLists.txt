rdkit_library(MolCatalog MolCatalogEntry.cpp MolCatalogParams.cpp
              LINK_LIBRARIES GraphMol Catalogs)

rdkit_test(testMolCatalog test1.cpp LINK_LIBRARIES MolCatalog FileParsers)

add_subdirectory(Wrap)



