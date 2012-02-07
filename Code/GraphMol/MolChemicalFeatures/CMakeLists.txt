rdkit_library(MolChemicalFeatures
              MolChemicalFeature.cpp MolChemicalFeatureDef.cpp 
              MolChemicalFeatureFactory.cpp FeatureParser.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse)

rdkit_headers(FeatureParser.h MolChemicalFeatureFactory.h
              MolChemicalFeatureDef.h  MolChemicalFeature.h 
              DEST GraphMol/MolChemicalFeatures)


rdkit_test(testFeatures testFeatures.cpp 
           LINK_LIBRARIES MolChemicalFeatures Catalogs
FileParsers SmilesParse SubstructMatch DataStructs
GraphMol RDGeometryLib RDGeneral
)

add_subdirectory(Wrap)

