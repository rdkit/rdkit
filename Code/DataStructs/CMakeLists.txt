if(RDK_OPTIMIZE_POPCNT)
  ADD_DEFINITIONS("-DRDK_OPTIMIZE_POPCNT")
endif()


rdkit_library(DataStructs
              BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp
              base64.cpp BitOps.cpp DiscreteDistMat.cpp
              DiscreteValueVect.cpp FPBReader.cpp MultiFPBReader.cpp
              LINK_LIBRARIES RDGeneral)
target_compile_definitions(DataStructs PRIVATE RDKIT_DATASTRUCTS_BUILD)

rdkit_headers(base64.h
              BitOps.h
              BitVect.h
              BitVects.h
              BitVectUtils.h
              DatastructsException.h
              DatastructsStreamOps.h
              DiscreteDistMat.h
              DiscreteValueVect.h
              ExplicitBitVect.h
              SparseBitVect.h
              SparseIntVect.h
              FPBReader.h
              MultiFPBReader.h
              DEST DataStructs)

rdkit_test(testDataStructs testDatastructs.cpp
           LINK_LIBRARIES DataStructs )

rdkit_test(testFPB testFPB.cpp
           LINK_LIBRARIES DataStructs )

rdkit_test(testMultiFPB testMultiFPB.cpp
           LINK_LIBRARIES DataStructs )

rdkit_catch_test(catchDataStructs catch_tests.cpp
           LINK_LIBRARIES DataStructs )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
