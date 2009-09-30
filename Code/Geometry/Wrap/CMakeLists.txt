rdkit_python_extension(rdGeometry
  Point.cpp rdGeometry.cpp
  DEST Geometry
  LINK_LIBRARIES RDGeometryLib RDGeneral RDBoost )
add_test( testPyGeometry python ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py )

