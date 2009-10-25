rdkit_library(MolCatalog MolCatalogEntry.cpp MolCatalogParams.cpp
              LINK_LIBRARIES GraphMol Catalogs)

rdkit_test(testMolCatalog test1.cpp LINK_LIBRARIES 
MolCatalog FileParsers Catalogs SmilesParse GraphMol DataStructs RDGeneral RDGeometryLib )

add_subdirectory(Wrap)



