if(NOT RDK_BUILD_YAEHMOP_SUPPORT)
  return()
endif(NOT RDK_BUILD_YAEHMOP_SUPPORT)

add_definitions(-DRDK_BUILD_YAEHMOP_SUPPORT)

if(NOT DEFINED YAEHMOP_DIR)
  set(YAEHMOP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/yaehmop")
endif()

if(NOT EXISTS "${YAEHMOP_DIR}/tightbind/bind.h")
    set(RELEASE_NO "2024.03.1")
    set(MD5 "ebbddca4f79ab71544cb1fef9a7eca8a")
    downloadAndCheckMD5("https://github.com/greglandrum/yaehmop/archive/refs/tags/v${RELEASE_NO}.tar.gz"
          "${CMAKE_CURRENT_SOURCE_DIR}/yaehmop-${RELEASE_NO}.tar.gz" ${MD5})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_CURRENT_SOURCE_DIR}/yaehmop-${RELEASE_NO}.tar.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    file(RENAME "yaehmop-${RELEASE_NO}" "${YAEHMOP_DIR}")
else()
  message("-- Found YAeHMOP source in ${YAEHMOP_DIR}")
endif()

set(yaehmop_INCLUDE_DIRS ${YAEHMOP_DIR}/..
     CACHE STRING "yaehmop Include File" FORCE)
include_directories(${yaehmop_INCLUDE_DIRS})

if(CMAKE_COMPILER_IS_GNUCXX AND NOT CYGWIN)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

# bring in the eHT code, but skip the targets (we'll get yaehmop_eht anyway because it's a dependency)
add_subdirectory(yaehmop/tightbind EXCLUDE_FROM_ALL True)

# set install dir for the yaehmop library:
INSTALL(TARGETS yaehmop_eht EXPORT rdkit-targets
        DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
        COMPONENT runtime )




set(EHT_PARAM_FILE ${CMAKE_CURRENT_SOURCE_DIR}/yaehmop/tightbind/eht_parms.dat )
install(FILES ${EHT_PARAM_FILE}
        DESTINATION ${RDKit_ShareDir}/Data
        COMPONENT data)

rdkit_library(EHTLib EHTTools.cpp SHARED LINK_LIBRARIES yaehmop_eht GraphMol )
target_compile_definitions(EHTLib PRIVATE RDKIT_EHTLIB_BUILD)
rdkit_headers(EHTTools.h DEST GraphMol)
rdkit_catch_test(testEHTLib1 test1.cpp
           LINK_LIBRARIES EHTLib FileParsers SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
