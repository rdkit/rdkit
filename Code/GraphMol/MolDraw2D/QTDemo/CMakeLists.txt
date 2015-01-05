
## for convenience
#find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost 1.46.0 COMPONENTS program_options iostreams filesystem )

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(CMAKE_COMPILER_IS_GNUCXX)
  ADD_DEFINITIONS("-fPIC")
endif()

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${Qt5Widgets_INCLUDE_DIRS}
  ${Qt5OpenGL_INCLUDE_DIRS}
  ${RDKIT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR} )

set( RDKitSV_SRCS ${RDKitSV_SOURCE_DIR}/QTDemo/rdkitsv_main.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/MolDisplay2DWidget.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVMainWindow.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVPanel.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVSettings.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/QTGet2Strings.cc
  ${RDKitSV_SOURCE_DIR}/QTDemo/QT4SelectItems.cc )

set( RDKitSV_INCS ${RDKitSV_SOURCE_DIR}/QTDemo/MolDisplay2DWidget.H
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVMainWindow.H
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVPanel.H
  ${RDKitSV_SOURCE_DIR}/QTDemo/RDKitSVSettings.H
  ${RDKitSV_SOURCE_DIR}/MolDraw2D.H
  ${RDKitSV_SOURCE_DIR}/MolDraw2DQt.H
  ${RDKitSV_SOURCE_DIR}/QTDemo/QTGet2Strings.H
  ${RDKitSV_SOURCE_DIR}/QTDemo/QT4SelectItems.H )

set( LIBS ${Qt5Widgets_LIBRARIES}
  ${Qt5OpenGL_LIBRARIES}
  ${Boost_LIBRARIES}
  ${RDKIT_LIBRARIES} z )

add_executable( rdkitsv ${RDKitSV_SRCS}
  ${RDKitSV_MOC_SRCS} ${RDKitSV_INCS} )
target_link_libraries( rdkitsv MolDraw2D ${LIBS}  )

