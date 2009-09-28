rdkit_python_extension(rdGeometry
                       Point.cpp UniformGrid3D.cpp rdGeometry.cpp
                       DEST Geometry
                       LINK_LIBRARIES RDGeometry RDBoost)

add_test(pyGeometry python ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py)

