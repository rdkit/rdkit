cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

option(RDK_BUILD_QT_SUPPORT "build support for QT drawing" ON )
option(RDK_BUILD_QT_DEMO "build the QT drawing demo" ON )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )

find_package( Boost 1.46.0)

# RDKit using my FindMyRDKit.cmake file in build_utils
# gcc 4.7.7 is really picky about the order of these libraries
find_package( MyRDKit COMPONENTS FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol RDGeometryLib )

if(CMAKE_COMPILER_IS_GNUCXX)
  ADD_DEFINITIONS("-fPIC")
endif()

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${RDKIT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR} )

set(EXECUTABLE_OUTPUT_PATH ${RDKitSV_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

if(RDK_BUILD_QT_SUPPORT)
  #find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
  find_package(Qt5Widgets REQUIRED)
  find_package(Qt5OpenGL REQUIRED)

include_directories(
  ${Qt5Widgets_INCLUDE_DIRS}
  ${Qt5OpenGL_INCLUDE_DIRS}
  )

  set(QTDRAW2D_SRC MolDraw2DQt.cc)
endif(RDK_BUILD_QT_SUPPORT)
  
add_library(MolDraw2D SHARED MolDraw2D.cc ${QTDRAW2D_SRC} )

if(RDK_BUILD_QT_DEMO)
  add_subdirectory(QTDemo)
endif(RDK_BUILD_QT_DEMO)
