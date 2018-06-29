rdkit_library(MolStandardize
              MolStandardize.cpp
	      Metal.cpp
	      Normalize.cpp
	      Validate.cpp
	      Charge.cpp
	      LINK_LIBRARIES RDGeneral GraphMol FragmentCatalog TransformCatalog
	      AcidBaseCatalog) 

rdkit_headers(MolStandardize.h
	      Metal.h
	      Normalize.h
	      Validate.h
	      Charge.h
	      DEST GraphMol/MolStandardize)
    
add_subdirectory(FragmentCatalog)
add_subdirectory(TransformCatalog)
add_subdirectory(AcidBaseCatalog)

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FragmentCatalog)
rdkit_test(molNormalizeTest testNormalize.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch TransformCatalog)
rdkit_test(molValidateTest testValidate.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers FragmentCatalog)
rdkit_test(molChargeTest testCharge.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers AcidBaseCatalog) 
