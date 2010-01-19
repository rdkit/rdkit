rdkit_python_extension(rdChemicalFeatures 
                       rdChemicalFeatures.cpp FreeChemicalFeature.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ChemicalFeatures GraphMol RDGeometryLib RDGeneral RDBoost) 

add_test(pyFeatures ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testFeatures.py)



