rdkit_library(MolStandardize
    Pipeline.cpp
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
        LINK_LIBRARIES ChemReactions ChemTransforms SmilesParse SubstructMatch Descriptors GraphMol )
target_compile_definitions(MolStandardize PRIVATE RDKIT_MOLSTANDARDIZE_BUILD)

rdkit_headers(
    Pipeline.h
    MolStandardize.h
    Metal.h
    Normalize.h
    Validate.h
    Charge.h
    Tautomer.h
    Fragment.h
        DEST GraphMol/MolStandardize)

rdkit_headers(
    FragmentCatalog/FragmentCatalogEntry.h
    FragmentCatalog/FragmentCatalogParams.h
    FragmentCatalog/FragmentCatalogUtils.h
        DEST GraphMol/MolStandardize/FragmentCatalog)

rdkit_headers(
    AcidBaseCatalog/AcidBaseCatalogEntry.h
    AcidBaseCatalog/AcidBaseCatalogParams.h
    AcidBaseCatalog/AcidBaseCatalogUtils.h
        DEST GraphMol/MolStandardize/AcidBaseCatalog)

rdkit_headers(
    TransformCatalog/TransformCatalogEntry.h
    TransformCatalog/TransformCatalogParams.h
    TransformCatalog/TransformCatalogUtils.h
        DEST GraphMol/MolStandardize/TransformCatalog)

rdkit_headers(
    TautomerCatalog/TautomerCatalogEntry.h
    TautomerCatalog/TautomerCatalogParams.h
    TautomerCatalog/TautomerCatalogUtils.h
        DEST GraphMol/MolStandardize/TautomerCatalog)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molNormalizeTest testNormalize.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molValidateTest testValidate.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molChargeTest testCharge.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molTautomerTest testTautomer.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molStandardizeSmallTest test2.cpp LINK_LIBRARIES MolStandardize )
rdkit_test(molFragmentTest testFragment.cpp LINK_LIBRARIES MolStandardize )
rdkit_catch_test(molStandardizeCatchTest catch_tests.cpp LINK_LIBRARIES MolStandardize )
rdkit_catch_test(molStandardizePipelineTest testPipeline.cpp LINK_LIBRARIES MolStandardize)


