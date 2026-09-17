remove_definitions(-DRDKIT_MOLDRAW2D_BUILD)

rdkit_nanobind_extension(rdMolDraw2D
                         rdMolDraw2D.cpp
                         DEST Chem/Draw
                         LINK_LIBRARIES MolDraw2D )

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolDraw2D.py)
