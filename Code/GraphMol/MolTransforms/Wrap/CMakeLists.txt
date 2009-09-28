rdkit_python_extension(rdMolTransforms 
                       rdMolTransforms.cpp
                       DEST Chem 
                       LINK_LIBRARIES MolTransforms RDBoost)

add_test(pyMolTransforms python ${CMAKE_CURRENT_SOURCE_DIR}/testMolTransforms.py)
