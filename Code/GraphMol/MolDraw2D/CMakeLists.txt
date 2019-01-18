
option(RDK_BUILD_QT_SUPPORT "build support for QT drawing" OFF )
option(RDK_BUILD_QT_DEMO "build the QT drawing demo" OFF )
option(RDK_BUILD_CAIRO_SUPPORT "build support for Cairo drawing" OFF )


rdkit_headers(MolDraw2D.h
  MolDraw2DSVG.h
  MolDraw2Dwx.h
  MolDraw2DUtils.h
  MolDraw2DDetails.h
  DEST GraphMol/MolDraw2D
)


remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLDRAW2D_BUILD)
rdkit_library(MolDraw2D MolDraw2D.cpp MolDraw2DSVG.cpp
  MolDraw2DDetails.cpp MolDraw2DUtils.cpp
  LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor MolTransforms RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol EigenSolvers
RDGeometryLib
  ${RDKit_THREAD_LIBS} )

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
  target_sources(MolDraw2D PRIVATE MolDraw2DCairo.cpp)
  rdkit_headers(MolDraw2DCairo.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_CAIRO_SUPPORT)

rdkit_test(moldraw2DTest1 test1.cpp LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol MolTransforms EigenSolvers
  RDGeometryLib
  MolDraw2D ${RDKit_THREAD_LIBS} )

rdkit_catch_test(moldraw2DTestCatch catch_tests.cpp LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol MolTransforms EigenSolvers
  RDGeometryLib
  MolDraw2D ${RDKit_THREAD_LIBS} )

rdkit_test(moldraw2DRxnTest1 rxn_test1.cpp LINK_LIBRARIES
    ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
    RDGeneral SubstructMatch Subgraphs GraphMol MolTransforms EigenSolvers
  RDGeometryLib
    MolDraw2D ${RDKit_THREAD_LIBS} ${EXTRA_LOCAL_LIBS} )

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)

add_subdirectory(Wrap)
