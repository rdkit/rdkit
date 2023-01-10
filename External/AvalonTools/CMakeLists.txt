if(NOT RDK_BUILD_AVALON_SUPPORT)
  return()
endif(NOT RDK_BUILD_AVALON_SUPPORT)

set(AVALON_VERSION "2.0.2")
if(NOT DEFINED AVALONTOOLS_DIR)
  set(AVALONTOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ava-formake-AvalonToolkit_${AVALON_VERSION}")
  set(fileToCheck "${AVALONTOOLS_DIR}/src/main/C/common/reaccsio.c")
  set(needDownload "TRUE")
  if(EXISTS "${fileToCheck}")
    set(needDownload "FALSE")
  endif()
else()
  string(REGEX REPLACE "\\\\" "/" AVALONTOOLS_DIR ${AVALONTOOLS_DIR})
  set(needDownload "FALSE")
endif()

set(AVALON_SRC_PATH "${AVALONTOOLS_DIR}/src/main/C")
set(AVALON_COMMON_PATH "${AVALON_SRC_PATH}/common")
set(AVALON_INCLUDE_PATH "${AVALON_SRC_PATH}/include")
set(AVALON_PROGRAMS_PATH "${AVALON_SRC_PATH}/programs")

if(needDownload)
  if(NOT DEFINED AVALONTOOLS_URL)
    set(AVALONTOOLS_URL "https://github.com/rohdebe1/ava-formake/archive/refs/tags/AvalonToolkit_${AVALON_VERSION}.tar.gz")
  endif()
  if(NOT DEFINED AVALONTOOLS_MD5SUM)
    set(AVALONTOOLS_MD5SUM "fc188383a8896802e948c977b73dbe71")
  endif()
  if(NOT DEFINED AVALONTOOLS_BASE)
    string(REGEX REPLACE "^.*/" "" AVALONTOOLS_BASE "${AVALONTOOLS_URL}")
  endif()
  downloadAndCheckMD5(${AVALONTOOLS_URL} "${CMAKE_CURRENT_SOURCE_DIR}/${AVALONTOOLS_BASE}" ${AVALONTOOLS_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf
    ${CMAKE_CURRENT_SOURCE_DIR}/AvalonToolkit_${AVALON_VERSION}.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set(reaccsio_c "${AVALONTOOLS_DIR}/src/main/C/common/reaccsio.c")
file(READ ${reaccsio_c} reaccsio_c_data)
if (NOT "${reaccsio_c_data}" MATCHES "fputc\\('\\\\0', fp\\)")
  string(REGEX REPLACE "( *)(int)( *)(Size;)" "\\1\\2\\3\\4\n\\1int\\3Len;" reaccsio_c_data "${reaccsio_c_data}")
  string(REGEX REPLACE "( *)(Size \\+= strlen\\(Line\\);)" "\\1Len = strlen(Line);\n\\1if (!Len) break;\n\\1Size += Len;" reaccsio_c_data "${reaccsio_c_data}")
  string(REGEX REPLACE "(bufsize = 5\\*80 \\+ mp->n_atoms\\*80 \\+ mp->n_bonds\\*80 \\+ mp->n_props\\*80)(;)" "\\1 + 1\\2" reaccsio_c_data "${reaccsio_c_data}")
  string(REGEX REPLACE "( *)(rewind\\(fp\\);)" "\\1fputc('\\\\0', fp);\n\\1fflush(fp);\n\\1\\2" reaccsio_c_data "${reaccsio_c_data}")
  file(WRITE ${reaccsio_c} "${reaccsio_c_data}")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-Wno-format-security -Wformat=0 -Wstringop-overflow=0 -Wformat-overflow=0 -Wno-unused-result )
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  add_compile_options(-Wno-absolute-value -Wno-return-type -Wno-tautological-overlap-compare)
endif()
if (MSVC)
  add_definitions("/D_CRT_SECURE_NO_WARNINGS")
  add_compile_options(/wd4018 /wd4101 /wd4146 /wd4334 /wd4477 /wd4715 /wd4716 /wd4996 /nologo)
endif(MSVC)

set(avalon_clib_srcs ${AVALON_COMMON_PATH}/layout.c
	${AVALON_COMMON_PATH}/symboltable.c
	${AVALON_COMMON_PATH}/patclean.c
	${AVALON_COMMON_PATH}/utilities.c
	${AVALON_COMMON_PATH}/symbol_lists.c
	${AVALON_COMMON_PATH}/stereo.c
	${AVALON_COMMON_PATH}/set.c
	${AVALON_COMMON_PATH}/perceive.c
	${AVALON_COMMON_PATH}/local.c
	${AVALON_COMMON_PATH}/graph.c
	${AVALON_COMMON_PATH}/geometry.c
	${AVALON_COMMON_PATH}/forio.c
	${AVALON_COMMON_PATH}/depictutil.c
	${AVALON_COMMON_PATH}/denormal.c
	${AVALON_COMMON_PATH}/casutils.c
	${AVALON_COMMON_PATH}/ssmatch.c
	${AVALON_COMMON_PATH}/rtutils.c
	${AVALON_COMMON_PATH}/smi2mol.c
	${AVALON_COMMON_PATH}/didepict.c
	${AVALON_COMMON_PATH}/pattern.c
	${AVALON_COMMON_PATH}/canonizer.c
	${AVALON_COMMON_PATH}/aacheck.c
	${AVALON_COMMON_PATH}/fixcharges.c
	${AVALON_PROGRAMS_PATH}/struchk.c
	${AVALON_COMMON_PATH}/reaccsio.c
	${AVALON_COMMON_PATH}/hashcode.c
  )

rdkit_library(avalon_clib ${avalon_clib_srcs})
target_compile_definitions(avalon_clib PRIVATE RDKIT_AVALONLIB_BUILD)
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(avalon_clib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${AVALON_INCLUDE_PATH})

rdkit_library(AvalonLib AvalonTools.cpp SHARED 
     LINK_LIBRARIES avalon_clib SubstructMatch FileParsers SmilesParse GraphMol DataStructs  )
target_compile_definitions(AvalonLib PRIVATE RDKIT_AVALONLIB_BUILD)
rdkit_headers(AvalonTools.h DEST GraphMol)
rdkit_test(testAvalonLib1 test1.cpp
           LINK_LIBRARIES AvalonLib SubstructMatch )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
