
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLCHEMICALFEATURES_BUILD)
rdkit_library(MolChemicalFeatures
              MolChemicalFeature.cpp MolChemicalFeatureDef.cpp 
              MolChemicalFeatureFactory.cpp FeatureParser.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse GraphMol)

rdkit_headers(FeatureParser.h MolChemicalFeatureFactory.h
              MolChemicalFeatureDef.h  MolChemicalFeature.h 
              DEST GraphMol/MolChemicalFeatures)


rdkit_test(testFeatures testFeatures.cpp 
           LINK_LIBRARIES MolChemicalFeatures FileParsers )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

