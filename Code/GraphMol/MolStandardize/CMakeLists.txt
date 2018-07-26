add_subdirectory(FragmentCatalog)
add_subdirectory(TransformCatalog)
add_subdirectory(AcidBaseCatalog)
add_subdirectory(TautomerCatalog)

rdkit_library(MolStandardize
              MolStandardize.cpp
	      Metal.cpp
	      Normalize.cpp
	      Validate.cpp
	      Charge.cpp
	      Tautomer.cpp
	      LINK_LIBRARIES RDGeneral GraphMol FragmentCatalog TransformCatalog
	      AcidBaseCatalog TautomerCatalog)

rdkit_headers(MolStandardize.h
	      Metal.h
	      Normalize.h
	      Validate.h
	      Charge.h
	      Tautomer.h
	      DEST GraphMol/MolStandardize)

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FragmentCatalog)
rdkit_test(molNormalizeTest testNormalize.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch TransformCatalog)
rdkit_test(molValidateTest testValidate.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers FragmentCatalog)
rdkit_test(molChargeTest testCharge.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers AcidBaseCatalog)
rdkit_test(molTautomerTest testTautomer.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	TautomerCatalog)
rdkit_test(molStandardizeSmallTest test2.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch)

find_package(Boost 1.56.0 COMPONENTS system iostreams REQUIRED)
if (NOT Boost_USE_STATIC_LIBS)
  add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
endif()

rdkit_test(compareMolVSTest testPCS.cpp LINK_LIBRARIES RDGeneral SmilesParse RDGeneral GraphMol MolStandardize  ${Boost_LIBRARIES})
