
PYTHON_ADD_MODULE(rdBase RDBase.cpp)
set_target_properties(rdBase PROPERTIES PREFIX "")
target_link_libraries(rdBase RDGeneral RDBoost ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
INSTALL(TARGETS rdBase
 LIBRARY DESTINATION ${RDKit_PythonDir} )
