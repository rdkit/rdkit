rdkit_library(DataStructs 
              BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp
              base64.cpp BitOps.cpp DiscreteDistMat.cpp DiscreteValueVect.cpp
              LINK_LIBRARIES RDGeneral)

rdkit_headers(base64.h
              BitOps.h
              BitVect.h
              BitVects.h
              BitVectUtils.h
              DatastructsException.h
              DiscreteDistMat.h
              DiscreteValueVect.h
              ExplicitBitVect.h
              SparseBitVect.h
              SparseIntVect.h DEST DataStructs)

rdkit_test(testDataStructs testDatastructs.cpp 
           LINK_LIBRARIES DataStructs RDGeneral)

add_subdirectory(Wrap)
