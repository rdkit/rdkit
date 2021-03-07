remove_definitions(-DRDKIT_MOLDRAW2DQT_BUILD)

if(RDK_BUILD_QT_SUPPORT)
  find_package(Qt5 COMPONENTS Core OpenGL REQUIRED)
  set( QtDrawingLibs Qt5::Core Qt5::OpenGL )
endif(RDK_BUILD_QT_SUPPORT)

rdkit_python_extension(rdMolDraw2DQt
                       rdMolDraw2DQt.cpp
                       DEST Chem/Draw
                       LINK_LIBRARIES MolDraw2DQt ${QtDrawingLibs} )

add_pytest(pyMolDraw2DQt ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2DQt.py)
