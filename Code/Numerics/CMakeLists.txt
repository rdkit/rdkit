add_subdirectory(Alignment)  
add_subdirectory(EigenSolvers)
add_subdirectory(Optimizer)

rdkit_test(testMatrices testMatrices.cpp LINK_LIBRARIES RDGeneral)
