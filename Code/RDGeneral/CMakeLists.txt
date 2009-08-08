add_library(RDGeneral Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp )

add_executable(testDict testDict.cpp)
target_link_libraries(testDict RDGeneral)

add_test(testDict ${EXECUTABLE_OUTPUT_PATH}/testDict)
