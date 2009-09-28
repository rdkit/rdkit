rdkit_library(Optimizer 
              BFGSOpt.cpp LinearSearch.cpp
              LINK_LIBRARIES RDGeometry)

rdkit_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer)


