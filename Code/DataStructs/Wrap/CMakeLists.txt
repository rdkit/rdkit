rdkit_python_extension(cDataStructs 
                       DataStructs.cpp DiscreteValueVect.cpp SparseIntVect.cpp 
                       wrap_SparseBV.cpp wrap_ExplicitBV.cpp wrap_BitOps.cpp 
                       wrap_Utils.cpp 
                       DEST DataStructs
                       LINK_LIBRARIES
                       RDGeneral DataStructs RDBoost)

add_test(pyBV python ${CMAKE_CURRENT_SOURCE_DIR}/testBV.py)
add_test(pyDiscreteValueVect python 
         ${CMAKE_CURRENT_SOURCE_DIR}/testDiscreteValueVect.py)
add_test(pySparseIntVect python
         ${CMAKE_CURRENT_SOURCE_DIR}/testSparseIntVect.py)
 



