rdkit_library(ChemicalFeatures FreeChemicalFeature.cpp)

rdkit_headers(ChemicalFeature.h
              FreeChemicalFeature.h DEST ChemicalFeatures)

rdkit_test(testChemicalFeatures testChemicalFeatures.cpp 
           LINK_LIBRARIES ChemicalFeatures RDGeneral)

add_subdirectory(Wrap)



