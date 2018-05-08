
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SUBGRAPHS_BUILD)
rdkit_library(Subgraphs Subgraphs.cpp SubgraphUtils.cpp
              LINK_LIBRARIES GraphMol )

rdkit_headers(Subgraphs.h
              SubgraphUtils.h DEST GraphMol/Subgraphs)

target_include_directories(Subgraphs PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings>)

rdkit_test(testSubgraphs1 test1.cpp LINK_LIBRARIES Subgraphs SmilesParse GraphMol RDGeometryLib RDGeneral  )
rdkit_test(testSubgraphs2 test2.cpp LINK_LIBRARIES SmilesParse Subgraphs GraphMol RDGeometryLib RDGeneral  )
