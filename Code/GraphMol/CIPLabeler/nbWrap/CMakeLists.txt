remove_definitions(-DRDKIT_CIPLABELER_BUILD)

rdkit_nanobind_extension(rdCIPLabeler rdCIPLabeler.cpp
                         DEST Chem
                         LINK_LIBRARIES GraphMol CIPLabeler)

add_pytest(pyCIPLabelsWrapTest ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/pyCIPLabelWrapTest.py)

add_pytest(pyCIPLabelsValidation ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/pyCIPLabelsValidation.py)
