if(CMAKE_COMPILER_IS_GNUCXX)
# Dict.cpp behaves badly if compiled with optimization under g++
SET_SOURCE_FILES_PROPERTIES(Dict.cpp PROPERTIES COMPILE_FLAGS -O0)
endif(CMAKE_COMPILER_IS_GNUCXX)

rdkit_library(RDGeneral 
              Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp)

rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

