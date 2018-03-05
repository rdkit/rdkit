rdkit_library(Optimizer 
              BFGSOpt.cpp LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib)

rdkit_headers(BFGSOpt.h DEST Numerics/Optimizer)

rdkit_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer RDGeometryLib RDGeneral )


