if(RDK_BUILD_FREESASA_SUPPORT)

if(NOT DEFINED FREESASA_DIR)
  set(FREESASA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/freesasa")
  set(needDownload "TRUE")
  if(EXISTS "${FREESASA_DIR}/src/freesasa.h")
    set(needDownload "FALSE")
  endif()
endif()

if(needDownload)
  # don't actually use the md5 here
  downloadAndCheckMD5("https://github.com/mittinatten/freesasa/archive/master.tar.gz"
        "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz" 
        "")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
    ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
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
    freesasa-master/src/classifier.c
    freesasa-master/src/classifier_naccess.c
    freesasa-master/src/classifier_oons.c
    freesasa-master/src/classifier_protor.c
    freesasa-master/src/coord.c
    freesasa-master/src/freesasa.c
    freesasa-master/src/lexer.c
    freesasa-master/src/nb.c
    freesasa-master/src/parser.c
    freesasa-master/src/pdb.c
    freesasa-master/src/rsa.c
    freesasa-master/src/sasa_lr.c
    freesasa-master/src/sasa_sr.c
    freesasa-master/src/selection.c
    freesasa-master/src/structure.c
    freesasa-master/src/util.c
    )

if(MSVC)
 SET_SOURCE_FILES_PROPERTIES(
    freesasa-master/src/classifier.c
    freesasa-master/src/classifier_naccess.c
    freesasa-master/src/classifier_oons.c
    freesasa-master/src/classifier_protor.c
    freesasa-master/src/coord.c
    freesasa-master/src/freesasa.c
    freesasa-master/src/lexer.c
    freesasa-master/src/nb.c
    freesasa-master/src/parser.c
    freesasa-master/src/pdb.c
    freesasa-master/src/rsa.c
    freesasa-master/src/sasa_lr.c
    freesasa-master/src/sasa_sr.c
    freesasa-master/src/selection.c
    freesasa-master/src/structure.c
    freesasa-master/src/util.c 
    PROPERTIES LANGUAGE CXX )
endif()

add_definitions(-DUSE_THREADS=0)
rdkit_library(freesasa_clib ${freesasa_clib_srcs})
include_directories("freesasa-master/src")

rdkit_library(FreeSASALib RDFreeSASA.cpp SHARED LINK_LIBRARIES freesasa_clib GraphMol DataStructs RDGeometryLib RDGeneral )
rdkit_headers(RDFreeSASA.h DEST GraphMol)

rdkit_test(testFreeSASALib testFreeSASA.cpp
           LINK_LIBRARIES FreeSASALib freesasa_clib FileParsers SmilesParse
                          GraphMol DataStructs RDGeometryLib RDGeneral)

endif(RDK_BUILD_FREESASA_SUPPORT)
    
