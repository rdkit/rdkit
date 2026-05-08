
rdkit_library(ReducedGraphs
              ReducedGraphs.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse GraphMol RDGeneral)
target_compile_definitions(ReducedGraphs PRIVATE RDKIT_REDUCEDGRAPHS_BUILD)

rdkit_headers(ReducedGraphs.h
              DEST GraphMol/ReducedGraphs)

rdkit_test(testReducedGraphs test1.cpp LINK_LIBRARIES
           ReducedGraphs )

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
add_subdirectory(boostWrap)
endif()

