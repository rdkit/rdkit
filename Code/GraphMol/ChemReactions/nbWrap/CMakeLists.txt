remove_definitions(-DRDKIT_CHEMREACTIONS_BUILD)
rdkit_nanobind_extension(rdChemReactions
                         Enumerate.cpp
                         rdChemReactions.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         MarvinParser
                         ChemReactions)

add_pytest(pyChemReactions
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testReactionWrapper.py)

add_pytest(pyChemReactionEnumerations
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testEnumerations.py)

add_pytest(pyChemReactionSanitize
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testSanitize.py)
