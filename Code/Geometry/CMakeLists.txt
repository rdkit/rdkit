
rdkit_library(RDGeometryLib 
              point.cpp Transform2D.cpp Transform3D.cpp 
              UniformGrid3D.cpp GridUtils.cpp
              LINK_LIBRARIES DataStructs RDGeneral)
target_compile_definitions(RDGeometryLib PRIVATE RDKIT_RDGEOMETRYLIB_BUILD)

rdkit_headers(Grid3D.h
              GridUtils.h
              point.h
              Transform2D.h
              Transform3D.h
              Transform.h
              UniformGrid3D.h
              Utils.h DEST Geometry)

rdkit_test(testTransforms testTransforms.cpp LINK_LIBRARIES RDGeometryLib )
rdkit_test(testGrid testGrid.cpp LINK_LIBRARIES RDGeometryLib)
rdkit_catch_test(geometryTestsCatch catch_tests.cpp
           LINK_LIBRARIES RDGeometryLib)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()