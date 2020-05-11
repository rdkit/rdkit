
option(RDK_BUILD_QT_SUPPORT "build support for QT drawing" OFF )
option(RDK_BUILD_QT_DEMO "build the QT drawing demo" OFF )
option(RDK_BUILD_CAIRO_SUPPORT "build support for Cairo drawing" OFF )
option(RDK_BUILD_FREETYPE_SUPPORT "build support for FreeType font handling" ON)

rdkit_headers(MolDraw2D.h
  MolDraw2DSVG.h
  MolDraw2Dwx.h
  MolDraw2DUtils.h
  DrawText.h
  DrawTextFT.h
  DrawTextCairo.h
  DrawTextSVG.h
  DEST GraphMol/MolDraw2D
)

remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLDRAW2D_BUILD)
rdkit_library(MolDraw2D MolDraw2D.cpp MolDraw2DSVG.cpp
  MolDraw2DDetails.cpp MolDraw2DUtils.cpp
  LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor MolTransforms 
  SubstructMatch Subgraphs GraphMol EigenSolvers )

if(RDK_BUILD_QT_SUPPORT)
  find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
  target_sources(MolDraw2D PRIVATE MolDraw2DQt.cpp)
  target_compile_options(MolDraw2D PRIVATE "-std=c++11")
  target_link_libraries(MolDraw2D PUBLIC Qt5::Widgets Qt5::OpenGL)
  rdkit_headers(MolDraw2DQt.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_QT_SUPPORT)

if(RDK_BUILD_CAIRO_SUPPORT)
  find_package(Cairo REQUIRED)
  target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
  target_link_libraries(MolDraw2D PUBLIC Cairo::Cairo)
  target_sources(MolDraw2D PRIVATE MolDraw2DCairo.cpp DrawTextCairo.cpp)
  rdkit_headers(MolDraw2DCairo.h DrawTextCairo.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_CAIRO_SUPPORT)

if(RDK_BUILD_FREETYPE_SUPPORT)
  find_package(Cairo REQUIRED)
  target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
  target_link_libraries(MolDraw2D PUBLIC Cairo::Cairo)
  find_package(Freetype REQUIRED)
  target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
  target_link_libraries(MolDraw2D PUBLIC Freetype::Freetype)
  target_sources(MolDraw2D PRIVATE MolDraw2DSVG.cpp MolDraw2DCairo.cpp DrawText.cpp
          DrawTextFT.cpp DrawTextCairo.cpp DrawTextSVG.cpp DrawTextFTCairo.cpp
          DrawTextFTSVG.cpp)
  rdkit_headers(MolDraw2DCairo.h MolDraw2DSVG.h DrawText.h DrawTextFT.h
          DrawTextCairo.h DrawTextSVG.h DrawTextFTCairo.h DrawTextFTSVG.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_FREETYPE_SUPPORT)

rdkit_test(moldraw2DTest1 test1.cpp LINK_LIBRARIES
  MolDraw2D )

rdkit_catch_test(moldraw2DTestCatch catch_tests.cpp LINK_LIBRARIES
  MolDraw2D )

rdkit_test(moldraw2DRxnTest1 rxn_test1.cpp LINK_LIBRARIES
  MolDraw2D )

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
