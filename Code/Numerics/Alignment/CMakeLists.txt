rdkit_library(Alignment AlignPoints.cpp LINK_LIBRARIES RDGeometryLib EigenSolvers)

rdkit_headers(AlignPoints.h DEST Numerics/Alignment)

rdkit_test(testAlignment testAlignment.cpp LINK_LIBRARIES Alignment EigenSolvers RDGeometryLib RDGeneral )

add_subdirectory(Wrap)
