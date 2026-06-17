remove_definitions(-DRDKIT_DESCRIPTORS_BUILD)
rdkit_nanobind_extension(rdMolDescriptors rdMolDescriptors.cpp
                         DEST Chem
                         LINK_LIBRARIES Descriptors Fingerprints)

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolDescriptors.py)
add_pytest(pyMolDescriptors3D ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/test3D.py)
