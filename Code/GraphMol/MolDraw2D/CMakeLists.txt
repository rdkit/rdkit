
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

if(RDK_BUILD_QT_SUPPORT)
  #find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
  find_package(Qt5Widgets REQUIRED)
  find_package(Qt5OpenGL REQUIRED)

  include_directories(
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5OpenGL_INCLUDE_DIRS}
    )

  set(QTDRAW2D_SRC MolDraw2DQt.cpp)
  set(EXTRA_LOCAL_LIBS ${EXTRA_LOCAL_LIBS} ${Qt5Widgets_LIBRARIES}  ${Qt5OpenGL_LIBRARIES} )
  rdkit_headers(MolDraw2DQt.h DEST GraphMol/MolDraw2D)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(RDK_BUILD_QT_SUPPORT)

if(RDK_BUILD_CAIRO_SUPPORT)
  find_package(Cairo REQUIRED)
  include_directories( ${CAIRO_INCLUDE_DIRS} )
  set(EXTRA_LOCAL_LIBS ${EXTRA_LOCAL_LIBS} ${CAIRO_LIBRARIES} )
  ADD_DEFINITIONS("-DRDK_CAIRO_BUILD")
  set(CAIRODRAW2D_SRC MolDraw2DCairo.cpp)
  rdkit_headers(MolDraw2DCairo.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_CAIRO_SUPPORT)

rdkit_library(MolDraw2D MolDraw2D.cpp MolDraw2DSVG.cpp
  ${CAIRODRAW2D_SRC} ${QTDRAW2D_SRC}
  MolDraw2DDetails.cpp MolDraw2DUtils.cpp
  LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor MolTransforms RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol EigenSolvers
RDGeometryLib
  ${RDKit_THREAD_LIBS} ${EXTRA_LOCAL_LIBS} )

rdkit_test(moldraw2DTest1 test1.cpp LINK_LIBRARIES
  ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol MolTransforms EigenSolvers
RDGeometryLib
  MolDraw2D ${RDKit_THREAD_LIBS} ${EXTRA_LOCAL_LIBS} )

rdkit_test(moldraw2DRxnTest1 rxn_test1.cpp LINK_LIBRARIES
    ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
    RDGeneral SubstructMatch Subgraphs GraphMol MolTransforms EigenSolvers
  RDGeometryLib
    MolDraw2D ${RDKit_THREAD_LIBS} ${EXTRA_LOCAL_LIBS} )

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)

add_subdirectory(Wrap)
