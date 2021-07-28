add_custom_target(ringdecomposerlib_support ALL)

if(NOT RDK_USE_URF)
  return()
endif(NOT RDK_USE_URF)

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

#if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32 AND (NOT RDK_INSTALL_STATIC_LIBS)))
if(WIN32)
    set (ringdecomposerlib_h ${URFLIB_DIR}/RingDecomposerLib.h)
    file(READ ${ringdecomposerlib_h} ringdecomposerlib_h_data)
    if (NOT "${ringdecomposerlib_h_data}" MATCHES "RDKIT_URFLIB_BUILD")
        string(REGEX REPLACE "(#if[ ]+\\([ ]*defined[ ]*\\([ ]*_WIN32[ ]*\\)[ ]*&&[ ]*)(defined\\([ ]*_MSC_VER[ ]*\\))" "\\1!defined( RDL_WIN_STATIC )" ringdecomposerlib_h_data "${ringdecomposerlib_h_data}")
        string(REGEX REPLACE "([ ]*)(#define[ ]+RDL_API[ ]+__declspec)(\\([ ]*dllexport[ ]*\\))" "\\1#ifdef RDKIT_URFLIB_BUILD\n\\1\\1\\2\\3\n\\1#else\n\\1\\1\\2(dllimport)\n\\1#endif" ringdecomposerlib_h_data "${ringdecomposerlib_h_data}")
        file(WRITE ${ringdecomposerlib_h} "${ringdecomposerlib_h_data}")
    endif()
endif()
rdkit_library(RingDecomposerLib ${URFSOURCES} SHARED)
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32 AND (NOT RDK_INSTALL_STATIC_LIBS)))
    target_compile_definitions(RingDecomposerLib PRIVATE RDKIT_URFLIB_BUILD)
endif()
install(TARGETS RingDecomposerLib DESTINATION ${RDKit_LibDir})
rdkit_headers(${URFLIB_DIR}/RingDecomposerLib.h DEST "")


set(RDK_URF_LIBS RingDecomposerLib
    CACHE STRING "the libraries for the URF calculation" FORCE)
