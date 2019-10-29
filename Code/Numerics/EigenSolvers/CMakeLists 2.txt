
add_definitions(-DRDKIT_EIGENSOLVERS_BUILD)
rdkit_library(EigenSolvers PowerEigenSolver.cpp LINK_LIBRARIES RDGeneral)

rdkit_headers(PowerEigenSolver.h DEST Numerics/EigenSolvers)

IF (LAPACK_FOUND)
rdkit_test(tesEigensolvers testEigenSolvers.cpp
           LINK_LIBRARIES EigenSolvers  RDGeneral )
target_include_directories(tesEigenSolvers PUBLIC ${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings)
ENDIF (LAPACK_FOUND)
