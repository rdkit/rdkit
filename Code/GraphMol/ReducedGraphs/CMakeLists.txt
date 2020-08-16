
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_REDUCEDGRAPHS_BUILD)
rdkit_library(ReducedGraphs
              ReducedGraphs.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse GraphMol RDGeneral)

rdkit_headers(ReducedGraphs.h
              DEST GraphMol/ReducedGraphs)

rdkit_test(testReducedGraphs test1.cpp LINK_LIBRARIES
           ReducedGraphs )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

