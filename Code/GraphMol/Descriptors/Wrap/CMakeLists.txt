rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       Descriptors Fingerprints RDBoost) 

add_test(pyMolDescriptors python ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)





