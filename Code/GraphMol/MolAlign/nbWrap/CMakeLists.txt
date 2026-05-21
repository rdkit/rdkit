remove_definitions(-DRDKIT_MOLALIGN_BUILD)
rdkit_nanobind_extension(rdMolAlign rdMolAlign.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         MolAlign O3AAlign DataStructs SmilesParse Descriptors)

add_pytest(pyMolAlign ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolAlign.py)





