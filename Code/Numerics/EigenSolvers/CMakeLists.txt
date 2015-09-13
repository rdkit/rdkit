rdkit_library(EigenSolvers PowerEigenSolver.cpp JacobiEigenSolver.cpp LINK_LIBRARIES RDGeneral)

rdkit_headers(PowerEigenSolver.h JacobiEigenSolver.h DEST Numerics/EigenSolvers)

IF (LAPACK_FOUND)
include_directories(${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings)
rdkit_test(testEigensolvers testEigenSolvers.cpp
           LINK_LIBRARIES EigenSolvers  RDGeneral )
ENDIF (LAPACK_FOUND)
