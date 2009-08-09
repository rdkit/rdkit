rdkit_python_extension(cDataStructs
  DataStructs.cpp DiscreteValueVect.cpp SparseIntVect.cpp 
  wrap_SparseBV.cpp wrap_ExplicitBV.cpp wrap_BitOps.cpp wrap_Utils.cpp
  DEST DataStructs
  LINK_LIBRARIES DataStructs RDGeneral RDBoost )

add_test( pyDataStructsTestBV python ${CMAKE_CURRENT_SOURCE_DIR}/testBV.py )
add_test( pyDataStructsTestDVV python ${CMAKE_CURRENT_SOURCE_DIR}/testDiscreteValueVect.py )
add_test( pyDataStructsTestSIV python ${CMAKE_CURRENT_SOURCE_DIR}/testSparseIntVect.py )


