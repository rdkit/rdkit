rdkit_python_extension(rdGeometry
  Point.cpp rdGeometry.cpp
  DEST Geometry
  LINK_LIBRARIES RDGeometry RDGeneral RDBoost )
add_test( pyGeometry python ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py )

