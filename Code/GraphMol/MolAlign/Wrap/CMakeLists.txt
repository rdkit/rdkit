rdkit_python_extension(rdMolAlign rdMolAlign.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       MolAlign RDBoost)

add_test(pyMolAlign python ${CMAKE_CURRENT_SOURCE_DIR}/testMolAlign.py)





