rdkit_python_extension(rdGeometry
                       Point.cpp UniformGrid3D.cpp rdGeometry.cpp
                       DEST Geometry
                       LINK_LIBRARIES RDGeometryLib DataStructs RDGeneral RDBoost)

add_pytest(testPyGeometry ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py)

