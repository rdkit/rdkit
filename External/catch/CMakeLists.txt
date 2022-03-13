add_custom_target(catch ALL)

if(NOT DEFINED CATCH_DIR)
  set(CATCH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/catch")
endif()
if(EXISTS "${CATCH_DIR}/catch.hpp")
  message("-- Found Catch2: ${CATCH_DIR}")
  set(CATCH_INCLUDE_DIR "${CATCH_DIR}" PARENT_SCOPE)
elseif(EXISTS "${CATCH_DIR}/single_include/catch2/catch.hpp")
  message("-- Found Catch2 source in ${CATCH_DIR}")
  set(CATCH_INCLUDE_DIR "${CATCH_DIR}/single_include/catch2" PARENT_SCOPE)
else()
  message("-- Catch not found in ${CATCH_DIR}")
  include(RDKitUtils)
  set(RELEASE_NO "2.13.8")
  downloadAndCheckMD5("https://github.com/catchorg/Catch2/archive/v${RELEASE_NO}.tar.gz"
        "${CMAKE_CURRENT_SOURCE_DIR}/v${RELEASE_NO}.tar.gz"
        "3639fb4cb8020de0dcf1fd5addb39418")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
    ${CMAKE_CURRENT_SOURCE_DIR}/v${RELEASE_NO}.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  file(RENAME "Catch2-${RELEASE_NO}" "${CATCH_DIR}")
  set(CATCH_INCLUDE_DIR "${CATCH_DIR}/single_include/catch2" PARENT_SCOPE)
endif()

