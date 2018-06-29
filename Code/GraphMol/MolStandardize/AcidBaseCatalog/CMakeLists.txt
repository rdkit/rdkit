rdkit_library(AcidBaseCatalog
	AcidBaseCatalogParams.cpp
	AcidBaseCatalogUtils.cpp
	AcidBaseCatalogEntry.cpp
	LINK_LIBRARIES Catalogs GraphMol RDGeometryLib RDGeneral Descriptors
	ChemTransforms)

rdkit_headers(AcidBaseCatalogParams.h
	AcidBaseCatalogUtils.h
	AcidBaseCatalogEntry.h
	DEST GraphMol/MolStandardize/AcidBaseCatalog)

#rdkit_test(acidBaseCatalogTest test1.cpp LINK_LIBRARIES AcidBaseCatalog
#	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
#	ChemTransforms Descriptors MolStandardize)
