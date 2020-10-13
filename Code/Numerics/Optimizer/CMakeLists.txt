
rdkit_library(Optimizer 
              BFGSOpt.cpp LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib Trajectory RDGeneral)
target_compile_definitions(Optimizer PRIVATE RDKIT_OPTIMIZER_BUILD)

rdkit_headers(BFGSOpt.h DEST Numerics/Optimizer)

rdkit_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer )


