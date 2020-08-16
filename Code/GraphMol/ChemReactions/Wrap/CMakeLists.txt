remove_definitions(-DRDKIT_CHEMREACTIONS_BUILD)
rdkit_python_extension(rdChemReactions
                       Enumerate.cpp
                       rdChemReactions.cpp
                       DEST Chem
                       LINK_LIBRARIES
ChemReactions )

add_pytest(pyChemReactions
         ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py)

add_pytest(pyChemReactionEnumerations
         ${CMAKE_CURRENT_SOURCE_DIR}/testEnumerations.py)

add_pytest(pyChemReactionSanitize
         ${CMAKE_CURRENT_SOURCE_DIR}/testSanitize.py)
