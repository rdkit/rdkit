rdkit_library(FragmentCatalog
	FragmentCalaogParams.cpp
	FragmentCalaogUtils.cpp
	FragmentCalaogEntry.cpp
	FragmentRemover.cpp
	LINK_LIBRARIES Catalogs GraphMol RDGeometryLib RDGeneral )

rdkit_headers(FragmentCalaogParams.h
	FragmentCalaogUtils.h
	FragmentCalaogEntry.h
	FragmentRemover.h	
	DEST GraphMol/MolStandardize/FragmentCatalog)

rdkit_test(fragmentCatalogTest test1.cpp LINK_LIBRARIES FragmentCatalog
	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	ChemTransforms)
