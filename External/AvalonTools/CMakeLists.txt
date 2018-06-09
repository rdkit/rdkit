if(RDK_BUILD_AVALON_SUPPORT)
if(NOT DEFINED AVALONTOOLS_DIR)
  set(AVALONTOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SourceDistribution")
  set(fileToPatch "${CMAKE_CURRENT_SOURCE_DIR}/SourceDistribution/common/reaccsio.c")
  set(needDownload "TRUE")
  if(EXISTS "${fileToPatch}")
    file(READ "${fileToPatch}" buffer)
    if("${buffer}" MATCHES "//MyFree\\(\\(char \\*\\)tempdir\\);")
      set(needDownload "FALSE")
    endif()
  endif()
else()
  string(REGEX REPLACE "\\\\" "/" AVALONTOOLS_DIR ${AVALONTOOLS_DIR})
  set(needDownload "FALSE")
endif()

set(AVALON_SRC_PATH ${AVALONTOOLS_DIR}/common)

if(needDownload)
  if(NOT DEFINED AVALONTOOLS_URL)
    set(AVALONTOOLS_URL "http://sourceforge.net/projects/avalontoolkit/files/AvalonToolkit_1.2/AvalonToolkit_1.2.0.source.tar")
  endif()
  if(NOT DEFINED AVALONTOOLS_MD5SUM)
    set(AVALONTOOLS_MD5SUM "092a94f421873f038aa67d4a6cc8cb54")
  endif()
  if(NOT DEFINED AVALONTOOLS_BASE)
    string(REGEX REPLACE "^.*/" "" AVALONTOOLS_BASE "${AVALONTOOLS_URL}")
  endif()
  downloadAndCheckMD5(${AVALONTOOLS_URL} "${CMAKE_CURRENT_SOURCE_DIR}/${AVALONTOOLS_BASE}" ${AVALONTOOLS_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf
    ${CMAKE_CURRENT_SOURCE_DIR}/AvalonToolkit_1.2.0.source.tar
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  # apply patch to AvalonTools
  configure_file("${fileToPatch}" "${fileToPatch}.orig" COPYONLY)
  file(READ "${fileToPatch}" buffer)
  string(REGEX REPLACE "MyFree\\(\\(char \\*\\)tempdir\\);"
    "//MyFree((char *)tempdir);" buffer "${buffer}")
  file(WRITE "${fileToPatch}" "${buffer}")
endif()

add_definitions(-DRDK_BUILD_AVALON_SUPPORT)
if (MSVC)
  add_definitions( "/D_CRT_SECURE_NO_WARNINGS /wd4224 /wd4101 /wd4018 /wd4996 /wd4244 /wd4305 /wd4013 /wd4146 /wd4334 /wd4715 /wd4715  /nologo" )
endif(MSVC)

set(avalon_clib_srcs ${AVALON_SRC_PATH}/layout.c
	${AVALON_SRC_PATH}/symboltable.c
	${AVALON_SRC_PATH}/patclean.c
	${AVALON_SRC_PATH}/utilities.c
	${AVALON_SRC_PATH}/symbol_lists.c
	${AVALON_SRC_PATH}/stereo.c
	${AVALON_SRC_PATH}/set.c
	${AVALON_SRC_PATH}/perceive.c
	${AVALON_SRC_PATH}/local.c
	${AVALON_SRC_PATH}/graph.c
	${AVALON_SRC_PATH}/geometry.c
	${AVALON_SRC_PATH}/forio.c
	${AVALON_SRC_PATH}/depictutil.c
	${AVALON_SRC_PATH}/denormal.c
	${AVALON_SRC_PATH}/casutils.c
	${AVALON_SRC_PATH}/ssmatch.c
	${AVALON_SRC_PATH}/rtutils.c
	${AVALON_SRC_PATH}/smi2mol.c
	${AVALON_SRC_PATH}/didepict.c
	${AVALON_SRC_PATH}/pattern.c
	${AVALON_SRC_PATH}/canonizer.c
	${AVALON_SRC_PATH}/aacheck.c
	${AVALON_SRC_PATH}/fixcharges.c
	${AVALON_SRC_PATH}/struchk.c
	${AVALON_SRC_PATH}/reaccsio.c
	${AVALON_SRC_PATH}/hashcode.c
  )

# we need this to ensure that builds continue
# to work on linux systems with older versions
# of glibc when we're building with gcc-4.1.
# Without this flag, we'll endup requiring
# glibc 2.7.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_definitions(-D_GNU_SOURCE=1)
endif()

     if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-parentheses -Wno-logical-op-parentheses -Wno-dangling-else -Wno-format")
     endif()


rdkit_library(avalon_clib ${avalon_clib_srcs})
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(avalon_clib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${AVALON_SRC_PATH})

add_definitions(-DRDKIT_AVALONLIB_BUILD)
rdkit_library(AvalonLib AvalonTools.cpp SHARED LINK_LIBRARIES avalon_clib FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral )
rdkit_headers(AvalonTools.h DEST GraphMol)
rdkit_test(testAvalonLib1 test1.cpp
           LINK_LIBRARIES AvalonLib avalon_clib
	   SubstructMatch
	   FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral)

add_subdirectory(Wrap)

endif(RDK_BUILD_AVALON_SUPPORT)
