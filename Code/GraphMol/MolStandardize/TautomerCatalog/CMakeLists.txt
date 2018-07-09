rdkit_library(TautomerCatalog
	TautomerCatalogParams.cpp
	TautomerCatalogUtils.cpp
	TautomerCatalogEntry.cpp
	LINK_LIBRARIES Catalogs GraphMol RDGeometryLib RDGeneral Descriptors
	ChemTransforms)

rdkit_headers(TautomerCatalogParams.h
	TautomerCatalogUtils.h
	TautomerCatalogEntry.h
	DEST GraphMol/MolStandardize/TautomerCatalog)

#rdkit_test(acidBaseCatalogTest test1.cpp LINK_LIBRARIES TautomerCatalog
#	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
#	ChemTransforms Descriptors MolStandardize)
