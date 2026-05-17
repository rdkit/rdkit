remove_definitions(-DRDKIT_MOLENUMERATOR_BUILD)
rdkit_nanobind_extension(rdMolEnumerator
                         rdMolEnumerator.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         MolEnumerator)

add_pytest(pyMolEnumerator ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
