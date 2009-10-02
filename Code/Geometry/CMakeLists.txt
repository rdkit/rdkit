rdkit_library(RDGeometryLib point.cpp Transform2D.cpp Transform3D.cpp 
              LINK_LIBRARIES RDGeneral)

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometryLib RDGeneral)

add_subdirectory(Wrap)
