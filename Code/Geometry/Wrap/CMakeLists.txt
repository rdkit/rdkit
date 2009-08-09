
PYTHON_ADD_MODULE(rdGeometry
  Point.cpp UniformGrid3D.cpp rdGeometry.cpp
)
set_target_properties(rdGeometry PROPERTIES PREFIX "")
target_link_libraries(rdGeometry 
  RDGeometry DataStructs RDGeneral RDBoost 
  ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
INSTALL(TARGETS rdGeometry
 LIBRARY DESTINATION ${RDKit_PythonDir}/Geometry )

