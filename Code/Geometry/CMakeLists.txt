rdkit_library(RDGeometryLib 
              point.cpp Transform2D.cpp Transform3D.cpp 
              UniformGrid3D.cpp GridUtils.cpp
              LINK_LIBRARIES DataStructs RDGeneral)

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometryLib)
rdkit_test(testGrid testGrid.cpp LINK_LIBRARIES RDGeometryLib)

add_subdirectory(Wrap)
