rdkit_library(TransformCatalog
	TransformCatalogParams.cpp
	TransformCatalogUtils.cpp
	TransformCatalogEntry.cpp
	#	TransformRemover.cpp
	LINK_LIBRARIES Catalogs GraphMol RDGeometryLib RDGeneral Descriptors
	ChemTransforms ChemReactions)

rdkit_headers(TransformCatalogParams.h
	TransformCatalogUtils.h
	TransformCatalogEntry.h
	#	TransformRemover.h	
	DEST GraphMol/MolStandardize/TransformCatalog)

#rdkit_test(transformCatalogTest test1.cpp LINK_LIBRARIES TransformCatalog
#	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
#	ChemTransforms Descriptors MolStandardize)
