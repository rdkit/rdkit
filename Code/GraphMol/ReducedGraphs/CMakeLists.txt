
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_REDUCEDGRAPHS_BUILD)
rdkit_library(ReducedGraphs
              ReducedGraphs.cpp
              LINK_LIBRARIES SubstructMatch SmilesParse GraphMol
                ${RDKit_THREAD_LIBS} )

rdkit_headers(ReducedGraphs.h
              DEST GraphMol/ReducedGraphs)

rdkit_test(testReducedGraphs test1.cpp LINK_LIBRARIES
           ReducedGraphs
           SubstructMatch SmilesParse
           GraphMol DataStructs RDGeometryLib
           RDGeneral  )

add_subdirectory(Wrap)

