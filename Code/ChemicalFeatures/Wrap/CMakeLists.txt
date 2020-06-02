remove_definitions(-DRDKIT_CHEMICALFEATURES_BUILD)
rdkit_python_extension(rdChemicalFeatures 
                       rdChemicalFeatures.cpp FreeChemicalFeature.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ChemicalFeatures GraphMol) 

add_pytest(pyFeatures ${CMAKE_CURRENT_SOURCE_DIR}/testFeatures.py)



