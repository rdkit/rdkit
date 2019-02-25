if(RDK_OPTIMIZE_NATIVE)
  ADD_DEFINITIONS("-DUSE_BUILTIN_POPCOUNT")
endif()


add_definitions(-DRDKIT_DATASTRUCTS_BUILD)
rdkit_library(DataStructs
              BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp
              base64.cpp BitOps.cpp DiscreteDistMat.cpp
              DiscreteValueVect.cpp FPBReader.cpp MultiFPBReader.cpp
              LINK_LIBRARIES RDGeneral)

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
           LINK_LIBRARIES DataStructs RDGeneral)

rdkit_test(testFPB testFPB.cpp
           LINK_LIBRARIES DataStructs RDGeneral)

rdkit_test(testMultiFPB testMultiFPB.cpp
           LINK_LIBRARIES DataStructs RDGeneral)

add_subdirectory(Wrap)
