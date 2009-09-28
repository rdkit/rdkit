include_directories(${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings)

rdkit_library(Subgraphs Subgraphs.cpp SubgraphUtils.cpp
              LINK_LIBRARIES GraphMol ${LAPACK_LIBRARIES})

rdkit_test(testSubgraphs1 test1.cpp LINK_LIBRARIES SmilesParse Subgraphs)
rdkit_test(testSubgraphs2 test2.cpp LINK_LIBRARIES SmilesParse Subgraphs)
