cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

set(RDBASE $ENV{RDBASE})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -std=c++11" )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system)

# specify where CMake can find the RDKit libraries
include_directories ( ${RDBASE}/Code )
link_directories ( ${RDBASE}/lib )

set(RDKit_LIBS RDKitChemReactions RDKitFileParsers RDKitSmilesParse RDKitDepictor
	       RDKitRDGeometryLib RDKitRDGeneral RDKitSubstructMatch RDKitSubgraphs
	       RDKitMolDraw2D RDKitGraphMol RDKitDistGeometry RDKitDistGeomHelpers
  	       RDKitMolAlign RDKitOptimizer RDKitForceField RDKitForceFieldHelpers
	       RDKitAlignment RDKitForceField  RDKitMolTransforms RDKitEigenSolvers )

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams ${RDKit_THREAD_LIBS} z  )

include_directories(${RDKIT_INCLUDE_DIR})
foreach( exnum RANGE 1 16 )
  message( "Example : ${exnum}" )
  set( exe "example${exnum}" )
  set( src "${exe}.cpp" )
  add_executable( ${exe} ${src} )
  target_link_libraries( ${exe} ${LIBS} ${RDKit_LIBS})
endforeach( exnum )

# add_executable( test_example test_example.cpp )
# target_link_libraries( test_example ${LIBS} )
