
rdkit_library(Subgraphs Subgraphs.cpp SubgraphUtils.cpp
              LINK_LIBRARIES GraphMol RDGeneral )
target_compile_definitions(Subgraphs PRIVATE RDKIT_SUBGRAPHS_BUILD)

rdkit_headers(Subgraphs.h
              SubgraphUtils.h DEST GraphMol/Subgraphs)

target_include_directories(Subgraphs PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings>)

rdkit_test(testSubgraphs1 test1.cpp LINK_LIBRARIES Subgraphs SmilesParse )
rdkit_test(testSubgraphs2 test2.cpp LINK_LIBRARIES SmilesParse Subgraphs )
