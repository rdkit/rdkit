remove_definitions(-DRDKIT_DESCRIPTORS_BUILD)

rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES Descriptors Fingerprints)

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)
add_pytest(pyMolDescriptors3D ${CMAKE_CURRENT_SOURCE_DIR}/test3D.py)
