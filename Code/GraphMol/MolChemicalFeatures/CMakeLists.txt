rdkit_library(MolChemicalFeatures
              MolChemicalFeature.cpp MolChemicalFeatureDef.cpp 
              MolChemicalFeatureFactory.cpp FeatureParser.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse)

rdkit_test(testFeatures testFeatures.cpp 
           LINK_LIBRARIES MolChemicalFeatures Catalogs
SmilesParse SubstructMatch DataStructs
GraphMol RDGeometryLib RDGeneral
)

add_subdirectory(Wrap)

