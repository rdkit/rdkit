rdkit_nanobind_extension(rdNormalizedDescriptors rdNormalizedDescriptors.cpp
                         DEST Chem
                         LINK_LIBRARIES NormalizedDescriptors Descriptors Fingerprints)

add_pytest(pyNormalizedDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testNormalizedDescriptors.py)

