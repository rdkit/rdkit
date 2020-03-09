cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system)

find_package( Cairo REQUIRED )

find_package( MyRDKit COMPONENTS ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs MolDraw2D GraphMol DistGeometry DistGeomHelpers
  MolAlign Optimizer ForceField ForceFieldHelpers Alignment
  ForceField  MolTransforms EigenSolvers )

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
  target_link_libraries( ${exe} ${LIBS} )
endforeach( exnum )

# add_executable( test_example test_example.cpp )
# target_link_libraries( test_example ${LIBS} )
