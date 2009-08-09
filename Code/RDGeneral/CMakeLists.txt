add_library(RDGeneral  SHARED Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )
INSTALL(TARGETS RDGeneral LIBRARY DESTINATION ${RDKit_LibDir})

rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

