cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )

## for convenience
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost 1.54.0 COMPONENTS program_options iostreams filesystem )

# RDKit using my FindMyRDKit.cmake file in build_utils
# gcc 4.7.7 is really picky about the order of these libraries
find_package( MyRDKit COMPONENTS FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol RDGeometryLib )

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${Qt5Widgets_INCLUDE_DIRS}
  ${Qt5OpenGL_INCLUDE_DIRS}
  ${RDKIT_INCLUDE_DIR} )

set(EXECUTABLE_OUTPUT_PATH ${RDKitSV_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

set( RDKitSV_SRCS ${RDKitSV_SOURCE_DIR}/rdkitsv_main.cc
  ${RDKitSV_SOURCE_DIR}/MolDisplay2DWidget.cc
  ${RDKitSV_SOURCE_DIR}/RDKitSVMainWindow.cc
  ${RDKitSV_SOURCE_DIR}/RDKitSVPanel.cc
  ${RDKitSV_SOURCE_DIR}/RDKitSVSettings.cc
  ${RDKitSV_SOURCE_DIR}/MolDraw2D.cc
  ${RDKitSV_SOURCE_DIR}/MolDraw2DQt.cc
  ${RDKitSV_SOURCE_DIR}/QTGet2Strings.cc
  ${RDKitSV_SOURCE_DIR}/QT4SelectItems.cc )

set( RDKitSV_INCS ${RDKitSV_SOURCE_DIR}/MolDisplay2DWidget.H
  ${RDKitSV_SOURCE_DIR}/RDKitSVMainWindow.H
  ${RDKitSV_SOURCE_DIR}/RDKitSVPanel.H
  ${RDKitSV_SOURCE_DIR}/RDKitSVSettings.H
  ${RDKitSV_SOURCE_DIR}/MolDraw2D.H
  ${RDKitSV_SOURCE_DIR}/MolDraw2DQt.H
  ${RDKitSV_SOURCE_DIR}/QTGet2Strings.H
  ${RDKitSV_SOURCE_DIR}/QT4SelectItems.H )

set( LIBS ${Qt5Widgets_LIBRARIES}
  ${Qt5OpenGL_LIBRARIES}
  ${Boost_LIBRARIES}
  ${RDKIT_LIBRARIES} z )

add_executable( rdkitsv ${RDKitSV_SRCS}
  ${RDKitSV_MOC_SRCS} ${RDKitSV_INCS} )
target_link_libraries( rdkitsv ${LIBS} )

