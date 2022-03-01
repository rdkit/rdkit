
rdkit_library(MolDraw2DQt MolDraw2DQt.cpp DrawTextQt.cpp
  LINK_LIBRARIES
  MolDraw2D)
target_compile_definitions(MolDraw2DQt PRIVATE RDKIT_MOLDRAW2DQT_BUILD)
target_compile_definitions(MolDraw2DQt PUBLIC "-DRDK_BUILD_QT_SUPPORT")
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
target_link_libraries(MolDraw2DQt PUBLIC Qt5::Widgets Qt5::OpenGL)

rdkit_headers(MolDraw2DQt.h DEST GraphMol/MolDraw2D/Qt)
if(RDK_INSTALL_STATIC_LIBS)
  target_link_libraries(MolDraw2DQt_static PUBLIC Qt5::Widgets Qt5::OpenGL)
  target_sources(MolDraw2DQt_static PRIVATE MolDraw2DQt.cpp DrawTextQt.cpp)
endif()


if(RDK_BUILD_FREETYPE_SUPPORT)
  target_compile_definitions(MolDraw2DQt PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
  target_sources(MolDraw2DQt PRIVATE DrawTextFTQt.cpp)
  if(RDK_INSTALL_STATIC_LIBS)
    target_compile_definitions(MolDraw2DQt_static PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
    target_sources(MolDraw2DQt_static PRIVATE DrawTextFTQt.cpp)
  endif()
  # I think Emscripten with Qt is nonsense, so 
  find_package(Freetype REQUIRED)
  target_link_libraries(MolDraw2DQt PUBLIC Freetype::Freetype)
  if(RDK_INSTALL_STATIC_LIBS)
    target_link_libraries(MolDraw2DQt_static PUBLIC Freetype::Freetype)
  endif()
endif(RDK_BUILD_FREETYPE_SUPPORT)

if(RDK_BUILD_QT_SUPPORT)
rdkit_catch_test(moldraw2DTestQt catch_qt.cpp LINK_LIBRARIES
  MolDraw2DQt Qt5::Widgets Qt5::OpenGL)
endif(RDK_BUILD_QT_SUPPORT)

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
