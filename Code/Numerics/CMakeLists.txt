add_subdirectory(Alignment)  
add_subdirectory(EigenSolvers)
add_subdirectory(Optimizer)
if(RDK_BUILD_ANI)
add_subdirectory(EigenSerializer)
endif(RDK_BUILD_ANI)
rdkit_headers(Matrix.h
              SquareMatrix.h
              SymmMatrix.h
              Vector.h 
              Conrec.h DEST Numerics)

rdkit_test(testMatrices testMatrices.cpp LINK_LIBRARIES RDGeneral)
rdkit_catch_test(testConrec testConrec.cpp LINK_LIBRARIES RDGeometryLib )
