
add_definitions(-DRDKIT_CHEMICALFEATURES_BUILD)
rdkit_library(ChemicalFeatures FreeChemicalFeature.cpp
              LINK_LIBRARIES RDGeneral RDGeometryLib)

rdkit_headers(ChemicalFeature.h
              FreeChemicalFeature.h DEST ChemicalFeatures)

rdkit_test(testChemicalFeatures testChemicalFeatures.cpp 
           LINK_LIBRARIES ChemicalFeatures RDGeneral)

add_subdirectory(Wrap)



