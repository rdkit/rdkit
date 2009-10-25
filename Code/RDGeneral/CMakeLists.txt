rdkit_library(RDGeneral 
              Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp)

rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

