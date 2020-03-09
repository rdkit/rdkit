
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLSTANDARDIZE_BUILD)
rdkit_library(MolStandardize
              MolStandardize.cpp
	      Metal.cpp
	      Normalize.cpp
	      Validate.cpp
	      Charge.cpp
	      Tautomer.cpp
	      Fragment.cpp
        FragmentCatalog/FragmentCatalogEntry.cpp
        FragmentCatalog/FragmentCatalogParams.cpp
        FragmentCatalog/FragmentCatalogUtils.cpp
        AcidBaseCatalog/AcidBaseCatalogEntry.cpp
        AcidBaseCatalog/AcidBaseCatalogParams.cpp
        AcidBaseCatalog/AcidBaseCatalogUtils.cpp
        TransformCatalog/TransformCatalogEntry.cpp
        TransformCatalog/TransformCatalogParams.cpp
        TransformCatalog/TransformCatalogUtils.cpp
        TautomerCatalog/TautomerCatalogEntry.cpp
        TautomerCatalog/TautomerCatalogParams.cpp
        TautomerCatalog/TautomerCatalogUtils.cpp
	      LINK_LIBRARIES ChemReactions ChemTransforms GraphMol RDGeneral )

rdkit_headers(MolStandardize.h
	      Metal.h
	      Normalize.h
	      Validate.h
	      Charge.h
	      Tautomer.h
	      Fragment.h
	      DEST GraphMol/MolStandardize)

rdkit_headers(FragmentCatalog/FragmentCatalogEntry.h
        FragmentCatalog/FragmentCatalogParams.h
        FragmentCatalog/FragmentCatalogUtils.h
	    DEST GraphMol/MolStandardize/FragmentCatalog)

rdkit_headers(AcidBaseCatalog/AcidBaseCatalogEntry.h
        AcidBaseCatalog/AcidBaseCatalogParams.h
        AcidBaseCatalog/AcidBaseCatalogUtils.h
	    DEST GraphMol/MolStandardize/AcidBaseCatalog)

rdkit_headers(TransformCatalog/TransformCatalogEntry.h
        TransformCatalog/TransformCatalogParams.h
        TransformCatalog/TransformCatalogUtils.h
	    DEST GraphMol/MolStandardize/TransformCatalog)

rdkit_headers(TautomerCatalog/TautomerCatalogEntry.h
        TautomerCatalog/TautomerCatalogParams.h
        TautomerCatalog/TautomerCatalogUtils.h
	    DEST GraphMol/MolStandardize/TautomerCatalog)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch )
rdkit_test(molNormalizeTest testNormalize.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch )
rdkit_test(molValidateTest testValidate.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers )
rdkit_test(molChargeTest testCharge.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers )
rdkit_test(molTautomerTest testTautomer.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch FileParsers )
rdkit_test(molStandardizeSmallTest test2.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch)
rdkit_test(molFragmentTest testFragment.cpp LINK_LIBRARIES MolStandardize
	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	ChemTransforms Descriptors)
rdkit_catch_test(molStandardizeCatchTest catch_tests.cpp LINK_LIBRARIES MolStandardize
	Catalogs SmilesParse RDGeneral GraphMol SubstructMatch FileParsers
	ChemTransforms Descriptors)

#find_package(Boost 1.56.0 COMPONENTS system iostreams REQUIRED)
#if (NOT Boost_USE_STATIC_LIBS)
#  add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
#endif()

# rdkit_test(compareMolVSTest testPCS.cpp LINK_LIBRARIES RDGeneral SmilesParse RDGeneral GraphMol MolStandardize  ${Boost_LIBRARIES})
