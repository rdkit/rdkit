
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS program_options iostreams filesystem)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(CMAKE_COMPILER_IS_GNUCXX)
  add_compile_options(-fPIC)
endif()

include_directories(
  ${RDKIT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR} )

set( RDKitSV_SRCS rdkitsv_main.cc
  MolDisplay2DWidget.cc
  RDKitSVMainWindow.cc
  RDKitSVPanel.cc
  RDKitSVSettings.cc
  QTGet2Strings.cc
  QT4SelectItems.cc )

set( RDKitSV_INCS MolDisplay2DWidget.H
  RDKitSVMainWindow.H
  RDKitSVPanel.H
  RDKitSVSettings.H
  QTGet2Strings.H
  QT4SelectItems.H )

set( LIBS Boost::program_options Boost::iostreams Boost::filesystem
  ${RDKit_THREAD_LIBS}
  FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol RDGeometryLib
   z QtDependencies)

add_executable( rdkitsv ${RDKitSV_SRCS}
  ${RDKitSV_MOC_SRCS} ${RDKitSV_INCS} )
target_link_libraries( rdkitsv MolDraw2D ${LIBS}  )
