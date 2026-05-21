remove_definitions(-DRDKIT_FORCEFIELDHELPERS_BUILD)
rdkit_nanobind_extension(rdForceFieldHelpers rdForceFields.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         ForceFieldHelpers)

add_pytest(pyForceFieldHelpers ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testHelpers.py)
