rdkit_library(EigenSolvers PowerEigenSolver.cpp LINK_LIBRARIES RDGeneral)

rdkit_headers(PowerEigenSolver.h DEST Numerics/EigenSolvers)

IF (LAPACK_FOUND)
include_directories(${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings)
rdkit_test(tesEigensolvers testEigenSolvers.cpp
           LINK_LIBRARIES EigenSolvers  RDGeneral )
ENDIF (LAPACK_FOUND)
