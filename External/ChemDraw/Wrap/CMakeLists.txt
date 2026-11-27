remove_definitions(-DRDKIT_CHEMDRAW_BUILD)
rdkit_python_extension(rdChemDraw
                       rdChemDraw.cpp
                       DEST Chem
                       LINK_LIBRARIES 
		       RDChemDrawLib RDChemDrawReactionLib
                      )

add_pytest(testChemDraw ${CMAKE_CURRENT_SOURCE_DIR}/testChemDraw.py)
