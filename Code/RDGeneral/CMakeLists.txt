if(MSVC)
add_library(dyRDGeneral SHARED Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )
INSTALL(TARGETS dyRDGeneral 
  LIBRARY DESTINATION ${RDKit_LibDir}
  ARCHIVE DESTINATION ${RDKit_LibDir}
  RUNTIME DESTINATION ${RDKit_LibDir}
)
add_library(RDGeneral Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )
INSTALL(TARGETS RDGeneral 
  LIBRARY DESTINATION ${RDKit_LibDir}
  ARCHIVE DESTINATION ${RDKit_LibDir}
  RUNTIME DESTINATION ${RDKit_LibDir}
)
else(MSVC)
add_library(RDGeneral Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )
INSTALL(TARGETS RDGeneral 
  LIBRARY DESTINATION ${RDKit_LibDir})
endif(MSVC)

rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

