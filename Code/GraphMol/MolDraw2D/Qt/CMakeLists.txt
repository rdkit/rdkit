
rdkit_library(MolDraw2DQt MolDraw2DQt.cpp DrawTextQt.cpp
  LINK_LIBRARIES
  MolDraw2D)
target_compile_definitions(MolDraw2DQt PRIVATE RDKIT_MOLDRAW2DQT_BUILD)
target_compile_definitions(MolDraw2DQt PUBLIC "-DRDK_BUILD_QT_SUPPORT")

add_library(QtDependencies INTERFACE)
if (RDK_USE_QT6)
  find_package(Qt6 COMPONENTS Core Widgets OpenGL REQUIRED)
  target_link_libraries(QtDependencies INTERFACE Qt6::Core Qt6::OpenGL Qt6::Widgets)
  set(RDK_QT_VERSION_STRING ${Qt6Core_VERSION})
else(RDK_USE_QT6)
  find_package(Qt5 COMPONENTS Core Widgets OpenGL REQUIRED)
  target_link_libraries(QtDependencies INTERFACE Qt5::Core Qt5::OpenGL Qt5::Widgets)
  set(RDK_QT_VERSION_STRING ${Qt5Core_VERSION})
  if (RDK_QT_VERSION_STRING STREQUAL "")
    # QT < 5.9 uses Qt5Core_VERSION_STRING
    set(RDK_QT_VERSION_STRING ${Qt5Core_VERSION_STRING})
  endif()
endif(RDK_USE_QT6)

message(STATUS "Building with Qt version ${RDK_QT_VERSION_STRING}")

target_link_libraries(MolDraw2DQt PUBLIC QtDependencies)
target_compile_definitions(MolDraw2DQt PRIVATE "-DRDK_QT_VERSION=\"${RDK_QT_VERSION_STRING}\"")

install(TARGETS QtDependencies EXPORT ${RDKit_EXPORTED_TARGETS} COMPONENT dev )

rdkit_headers(MolDraw2DQt.h DEST GraphMol/MolDraw2D/Qt)
if(RDK_INSTALL_STATIC_LIBS)
  target_link_libraries(MolDraw2DQt_static PUBLIC QtDependencies)
  target_compile_definitions(MolDraw2DQt_static PRIVATE "-DRDK_QT_VERSION=\"${RDK_QT_VERSION_STRING}\"")
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
  rdkit_catch_test(moldraw2DTestQt catch_qt.cpp LINK_LIBRARIES MolDraw2DQt)
endif(RDK_BUILD_QT_SUPPORT)

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
