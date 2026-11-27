rdkit_nanobind_extension(rdDASHTree
                         rdDASHTree.cpp
                         DEST Chem
                         LINK_LIBRARIES DASHTree)
add_pytest(pyDASHTree ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testDASHTree.py)
