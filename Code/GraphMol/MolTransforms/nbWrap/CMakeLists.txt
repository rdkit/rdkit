remove_definitions(-DRDKIT_MOLTRANSFORMS_BUILD)
rdkit_nanobind_extension(rdMolTransforms
                         rdMolTransforms.cpp
                         DEST Chem
                         LINK_LIBRARIES MolTransforms)

add_pytest(pyMolTransforms ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolTransforms.py)
