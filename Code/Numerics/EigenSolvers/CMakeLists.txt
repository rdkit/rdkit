rdkit_library(EigenSolvers PowerEigenSolver.cpp LINK_LIBRARIES RDGeneral)

IF (LAPACK_FOUND)
include_directories(${CMAKE_SOURCE_DIR}/External/boost-numeric-bindings)
rdkit_test(tesEigensolvers testEigenSolvers.cpp
           LINK_LIBRARIES EigenSolvers ${LAPACK_LIBRARIES} RDGeneral )
ENDIF (LAPACK_FOUND)
