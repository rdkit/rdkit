remove_definitions(-DRDKIT_MOLDRAW2D_BUILD)

if(RDK_BUILD_QT_SUPPORT)
  find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
  set( QtDrawingLibs Qt5::Widgets Qt5::OpenGL )
endif(RDK_BUILD_QT_SUPPORT)

rdkit_python_extension(rdMolDraw2D
                       rdMolDraw2D.cpp
                       DEST Chem/Draw
                       LINK_LIBRARIES MolDraw2D ${QtDrawingLibs} )

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2D.py)
if(RDK_BUILD_QT_SUPPORT)
  add_pytest(pyMolDraw2DQt ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2DQt.py)
endif(RDK_BUILD_QT_SUPPORT)
