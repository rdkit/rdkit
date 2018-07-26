
rdkit_test(fragmentCatalogTest test1.cpp LINK_LIBRARIES MolStandardize
	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	ChemTransforms Descriptors MolStandardize)
