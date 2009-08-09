add_library(RDGeometry point.cpp Transform2D.cpp Transform3D.cpp UniformGrid3D.cpp GridUtils.cpp )

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometry DataStructs RDGeneral)
rdkit_test(testGrid testGrid.cpp LINK_LIBRARIES RDGeometry DataStructs RDGeneral)

add_subdirectory(Wrap)
