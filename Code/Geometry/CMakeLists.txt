rdkit_library(RDGeometryLib 
              point.cpp Transform2D.cpp Transform3D.cpp 
              UniformGrid3D.cpp GridUtils.cpp
              LINK_LIBRARIES DataStructs RDGeneral)

rdkit_headers(Grid3D.h
              GridUtils.h
              point.h
              Transform2D.h
              Transform3D.h
              Transform.h
              UniformGrid3D.h
              Utils.h DEST Geometry)

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometryLib DataStructs RDGeneral )
rdkit_test(testGrid testGrid.cpp LINK_LIBRARIES RDGeometryLib DataStructs RDGeneral )

add_subdirectory(Wrap)
