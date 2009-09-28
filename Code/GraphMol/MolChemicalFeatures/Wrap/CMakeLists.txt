rdkit_python_extension(rdMolChemicalFeatures
                       rdMolChemicalFeatures.cpp MolChemicalFeature.cpp 
                       MolChemicalFeatureFactory.cpp ChemicalFeatureUtils.cpp
		       DEST Chem
		       LINK_LIBRARIES MolChemicalFeatures RDBoost)

add_test(pyChemicalFeatures python
         ${CMAKE_CURRENT_SOURCE_DIR}/testChemicalFeatures.py)




