remove_definitions(-DRDKIT_MOLINTERCHANGE_BUILD)
rdkit_nanobind_extension(rdMolInterchange
                         rdMolInterchange.cpp
                         DEST Chem
                         LINK_LIBRARIES MolInterchange)

add_pytest(pyMolInterchange ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolInterchange.py)
