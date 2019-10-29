
add_definitions(-DRDKIT_OPTIMIZER_BUILD)
rdkit_library(Optimizer 
              BFGSOpt.cpp LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib Trajectory)

rdkit_headers(BFGSOpt.h DEST Numerics/Optimizer)

rdkit_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer RDGeometryLib RDGeneral )


