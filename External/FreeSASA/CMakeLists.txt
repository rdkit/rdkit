if(RDK_BUILD_FREESASA_SUPPORT)

add_definitions("-DRDKIT_FREESASALIB_BUILD")

if(NOT DEFINED FREESASA_DIR)
  set(FREESASA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/freesasa")
  set(needDownload "TRUE")
  if(EXISTS "${FREESASA_DIR}/src/freesasa.h")
    set(needDownload "FALSE")
  endif()
endif()

if(needDownload)
  # don't actually use the md5 here
  set(FREESASA_VERSION "2.0.3")
  set(FREESASA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/freesasa-${FREESASA_VERSION}")
  if(NOT EXISTS "${FREESASA_SRC_DIR}/src")
    downloadAndCheckMD5("https://github.com/mittinatten/freesasa/releases/download/${FREESASA_VERSION}/freesasa-${FREESASA_VERSION}.tar.gz"
          "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
          "")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    message("Copying freesasa2.c to ${FREESASA_SRC_DIR}/src/freesasa2.c")
    file(COPY freesasa2.c DESTINATION ${FREESASA_SRC_DIR}/src/)
  else()
    message("FreeSASA found in ${FREESASA_SRC_DIR}, no download required")
  endif()
endif()

macro(use_c99)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")
    endif ()
  else ()
    set (CMAKE_C_STANDARD 99)
  endif ()
endmacro(use_c99)

use_c99()

set (freesasa_clib_srcs
    ${FREESASA_SRC_DIR}/src/classifier.c
    ${FREESASA_SRC_DIR}/src/classifier_naccess.c
    ${FREESASA_SRC_DIR}/src/classifier_oons.c
    ${FREESASA_SRC_DIR}/src/classifier_protor.c
    ${FREESASA_SRC_DIR}/src/coord.c
    ${FREESASA_SRC_DIR}/src/freesasa2.c
    ${FREESASA_SRC_DIR}/src/lexer.c
    ${FREESASA_SRC_DIR}/src/log.c
    ${FREESASA_SRC_DIR}/src/nb.c
    ${FREESASA_SRC_DIR}/src/node.c
    ${FREESASA_SRC_DIR}/src/parser.c
    ${FREESASA_SRC_DIR}/src/pdb.c
    ${FREESASA_SRC_DIR}/src/rsa.c
    ${FREESASA_SRC_DIR}/src/sasa_lr.c
    ${FREESASA_SRC_DIR}/src/sasa_sr.c
    ${FREESASA_SRC_DIR}/src/selection.c
    ${FREESASA_SRC_DIR}/src/structure.c
    ${FREESASA_SRC_DIR}/src/util.c
    )
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set(freesasa_additional_exports ${CMAKE_CURRENT_SOURCE_DIR}/additional_exports.def)
  file(WRITE ${freesasa_additional_exports}
    "EXPORTS\n"
    "freesasa_default_parameters DATA\n"
    "freesasa_protor_classifier DATA\n"
    "freesasa_naccess_classifier DATA\n"
    "freesasa_oons_classifier DATA\n"
  )
  set (freesasa_clib_srcs
    ${freesasa_clib_srcs}
    ${freesasa_additional_exports}
  )
endif()

set (freesasa_h ${FREESASA_SRC_DIR}/src/freesasa.h)
file(READ ${freesasa_h} freesasa_h_data)
string(REGEX REPLACE "(#include <stdio.h>)" "\\1\n#include <RDGeneral/export.h>" freesasa_h_data "${freesasa_h_data}")
string(REGEX REPLACE "([^R][^D][^K][^I][^T][^_][^F][^R][^E][^E][^S][^A][^S][^A][^_][^C][^L][^I][^B][^_][^E][^X][^P][^O][^R][^T][^ ])(extern const)" "\\1RDKIT_FREESASA_CLIB_EXPORT \\2" freesasa_h_data "${freesasa_h_data}")
file(WRITE ${freesasa_h} "${freesasa_h_data}")

add_definitions(-DUSE_THREADS=0)
add_definitions(-DUSE_JSON=0)
add_definitions(-DUSE_XML=0)
rdkit_library(freesasa_clib ${freesasa_clib_srcs})
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(freesasa_clib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()
include_directories("${FREESASA_SRC_DIR}/src")

rdkit_library(FreeSASALib RDFreeSASA.cpp SHARED LINK_LIBRARIES freesasa_clib GraphMol DataStructs RDGeometryLib RDGeneral )
rdkit_headers(RDFreeSASA.h DEST GraphMol)

rdkit_test(testFreeSASALib testFreeSASA.cpp
           LINK_LIBRARIES FreeSASALib freesasa_clib FileParsers SmilesParse
                          GraphMol DataStructs RDGeometryLib RDGeneral)

add_subdirectory(Wrap)

endif(RDK_BUILD_FREESASA_SUPPORT)
