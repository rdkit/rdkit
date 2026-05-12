remove_definitions(-DRDKIT_MOLENUMERATOR_BUILD)
rdkit_python_extension(rdMolEnumerator
                       rdMolEnumerator.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       MolEnumerator )

add_pytest(pyMolEnumerator ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
