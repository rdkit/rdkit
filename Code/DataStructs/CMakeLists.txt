rdkit_library(DataStructs 
              BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp
              base64.cpp BitOps.cpp DiscreteDistMat.cpp DiscreteValueVect.cpp
              LINK_LIBRARIES RDGeneral)

rdkit_test(testDataStructs testDatastructs.cpp 
           LINK_LIBRARIES DataStructs RDGeneral)

add_subdirectory(Wrap)
