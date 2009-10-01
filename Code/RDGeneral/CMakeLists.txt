add_library(RDGeneral Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )
INSTALL(TARGETS RDGeneral 
  DESTINATION ${RDKit_LibDir})
#endif(MSVC)

rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

