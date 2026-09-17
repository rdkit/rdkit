remove_definitions(-DRDKIT_MOLDRAW2DQT_BUILD)

rdkit_nanobind_extension(rdMolDraw2DQt
                         rdMolDraw2DQt.cpp
                         DEST Chem/Draw
                         LINK_LIBRARIES MolDraw2DQt)

add_pytest(pyMolDraw2DQt ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolDraw2DQt.py)
