rdkit_python_extension(rdForceFieldHelpers rdForceFields.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ForceFieldHelpers RDBoost) 

add_test(pyForceFieldHelpers python ${CMAKE_CURRENT_SOURCE_DIR}/testHelpers.py)




