rdkit_library(FragmentCatalog
	FragmentCatalogParams.cpp
	FragmentCatalogUtils.cpp
	FragmentCatalogEntry.cpp
	FragmentRemover.cpp
	LINK_LIBRARIES Catalogs GraphMol RDGeometryLib RDGeneral Descriptors
	ChemTransforms)

rdkit_headers(FragmentCatalogParams.h
	FragmentCatalogUtils.h
	FragmentCatalogEntry.h
	FragmentRemover.h	
	DEST GraphMol/MolStandardize/FragmentCatalog)

rdkit_test(fragmentCatalogTest test1.cpp LINK_LIBRARIES FragmentCatalog
	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	ChemTransforms Descriptors MolStandardize)
