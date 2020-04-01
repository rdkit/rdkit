remove_definitions(-DRDKIT_CIPLABELLING_BUILD)

rdkit_python_extension(rdCIPLabelling rdCIPLabelling.cpp
                       DEST Chem
                       LINK_LIBRARIES GraphMol CIPLabelling)

add_pytest(pyCIPLabelsValidation ${CMAKE_CURRENT_SOURCE_DIR}/pyCIPLabelsValidation.py)
