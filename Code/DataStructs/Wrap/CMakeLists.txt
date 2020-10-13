remove_definitions(-DRDKIT_DATASTRUCTS_BUILD)
rdkit_python_extension(cDataStructs
                       DataStructs.cpp DiscreteValueVect.cpp SparseIntVect.cpp
                       wrap_SparseBV.cpp wrap_ExplicitBV.cpp wrap_BitOps.cpp
                       wrap_FPB.cpp
                       wrap_Utils.cpp
                       DEST DataStructs
                       LINK_LIBRARIES
                       DataStructs)

add_pytest(pyBV
         ${CMAKE_CURRENT_SOURCE_DIR}/testBV.py)
add_pytest(pyDiscreteValueVect
         ${CMAKE_CURRENT_SOURCE_DIR}/testDiscreteValueVect.py)
add_pytest(pySparseIntVect
         ${CMAKE_CURRENT_SOURCE_DIR}/testSparseIntVect.py)
add_pytest(pyFPB
         ${CMAKE_CURRENT_SOURCE_DIR}/testFPB.py)
