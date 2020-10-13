
rdkit_library(ChemicalFeatures FreeChemicalFeature.cpp
              LINK_LIBRARIES RDGeneral RDGeometryLib)
target_compile_definitions(ChemicalFeatures PRIVATE RDKIT_CHEMICALFEATURES_BUILD)

rdkit_headers(ChemicalFeature.h
              FreeChemicalFeature.h DEST ChemicalFeatures)

rdkit_test(testChemicalFeatures testChemicalFeatures.cpp 
           LINK_LIBRARIES ChemicalFeatures)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()


