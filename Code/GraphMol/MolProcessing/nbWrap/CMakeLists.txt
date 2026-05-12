rdkit_nanobind_extension(rdMolProcessing
                         rdMolProcessing.cpp
                         DEST Chem
                         LINK_LIBRARIES MolProcessing GraphMol)

add_pytest(pyMolProcessing ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolProcessing.py)
