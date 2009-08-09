
PYTHON_ADD_MODULE(cDataStructs DataStructs.cpp DiscreteValueVect.cpp SparseIntVect.cpp 
  wrap_SparseBV.cpp wrap_ExplicitBV.cpp wrap_BitOps.cpp wrap_Utils.cpp )
set_target_properties(cDataStructs PROPERTIES PREFIX "")
target_link_libraries(cDataStructs DataStructs RDGeneral RDBoost ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
INSTALL(TARGETS cDataStructs 
 LIBRARY DESTINATION ${RDKit_PythonDir}/DataStructs )

