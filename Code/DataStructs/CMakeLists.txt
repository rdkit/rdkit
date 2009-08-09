add_library(DataStructs BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp 
	base64.cpp BitOps.cpp DiscreteDistMat.cpp DiscreteValueVect.cpp )

add_executable(testDatastructs testDatastructs.cpp)
target_link_libraries(testDatastructs DataStructs RDGeneral)
add_test(testDatastructs ${EXECUTABLE_OUTPUT_PATH}/testDatastructs)

add_subdirectory(Wrap)

