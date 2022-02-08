remove_definitions(-DRDKIT_MOLDRAW2DQT_BUILD)

find_package(Qt5 COMPONENTS Core OpenGL REQUIRED)

rdkit_python_extension(rdMolDraw2DQt
                       rdMolDraw2DQt.cpp
                       DEST Chem/Draw
                       LINK_LIBRARIES MolDraw2DQt Qt5::Core Qt5::OpenGL )

add_pytest(pyMolDraw2DQt ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2DQt.py)
