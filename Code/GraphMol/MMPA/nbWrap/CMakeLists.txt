remove_definitions(-DRDKIT_MMPA_BUILD)
rdkit_nanobind_extension(rdMMPA
                         rdMMPA.cpp
                         DEST Chem
                         LINK_LIBRARIES MMPA SmilesParse)

add_pytest(pyMMPA ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMMPA.py)
