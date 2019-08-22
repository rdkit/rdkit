add_custom_target(ringdecomposerlib_support ALL)

if(RDK_USE_URF)

if(NOT DEFINED URFLIB_DIR)
    set(URFLIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/RingDecomposerLib/src/RingDecomposerLib")
endif()

if(NOT EXISTS "${URFLIB_DIR}/RingDecomposerLib.h")
    set(RELEASE_NO "1.1.3_rdkit")
    set(MD5 "e9a0bcdda8b921a35e812b9888a9a874")
    downloadAndCheckMD5("https://github.com/rareylab/RingDecomposerLib/archive/v${RELEASE_NO}.tar.gz"
            "${CMAKE_CURRENT_SOURCE_DIR}/RingDecomposerLib-v${RELEASE_NO}.tar.gz" ${MD5})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
        ${CMAKE_CURRENT_SOURCE_DIR}/RingDecomposerLib-v${RELEASE_NO}.tar.gz
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(RENAME "RingDecomposerLib-${RELEASE_NO}" "RingDecomposerLib")
else()
    message("-- Found RingDecomposerLib source in ${URFLIB_DIR}")
endif()

set(urflib_INCLUDE_DIRS ${URFLIB_DIR}
    CACHE STRING "RingDecomposerLib Include Dir" FORCE)
file(GLOB URFSOURCES "${URFLIB_DIR}/*.c")

rdkit_library(RingDecomposerLib ${URFSOURCES} SHARED)
install(TARGETS RingDecomposerLib DESTINATION ${RDKit_LibDir})
rdkit_headers(${URFLIB_DIR}/RingDecomposerLib.h DEST "")


rdkit_library(URFLib URF.cpp SHARED LINK_LIBRARIES RingDecomposerLib )
set(RDK_URF_LIBS URFLib RingDecomposerLib
    CACHE STRING "the libraries for the URF calculation" FORCE)

endif(RDK_USE_URF)
