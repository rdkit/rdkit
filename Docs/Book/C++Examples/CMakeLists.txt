cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11" )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system thread)

find_package( MyRDKit COMPONENTS FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs GraphMol DistGeometry DistGeomHelpers
  MolAlign Optimizer RDGeometryLib ForceField ForceFieldHelpers Alignment
  ForceField ForceFieldHelpers MolTransforms EigenSolvers)

include_directories( ${RDKIT_INCLUDE_DIR} )

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

set( LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES} z )

add_executable( example1 example1.cpp )
target_link_libraries( example1 ${LIBS} )

add_executable( example2 example2.cpp )
target_link_libraries( example2 ${LIBS} )

add_executable( example3 example3.cpp )
target_link_libraries( example3 ${LIBS} )

add_executable( example4 example4.cpp )
target_link_libraries( example4 ${LIBS} )

add_executable( example5 example5.cpp )
target_link_libraries( example5 ${LIBS} )

add_executable( example6 example6.cpp )
target_link_libraries( example6 ${LIBS} )

add_executable( example7 example7.cpp )
target_link_libraries( example7 ${LIBS} )
