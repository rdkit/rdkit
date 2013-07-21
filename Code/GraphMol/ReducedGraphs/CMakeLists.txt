rdkit_library(ReducedGraphs
              ReducedGraphs.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse GraphMol
                ${RDKit_THREAD_LIBS} )

rdkit_headers(ReducedGraphs.h
              DEST GraphMol/ReducedGraphs)

rdkit_test(testReducedGraphs test1.cpp LINK_LIBRARIES
           ReducedGraphs_static
           SubstructMatch SmilesParse
           Subgraphs GraphMol DataStructs RDGeometryLib
           RDGeneral  )


