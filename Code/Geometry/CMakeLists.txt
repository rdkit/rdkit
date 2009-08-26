add_library(RDGeometry point.cpp Transform2D.cpp Transform3D.cpp )

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometry RDGeneral)

add_subdirectory(Wrap)
