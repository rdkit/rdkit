if(NOT RDK_BUILD_ROSHAMBOSHAPE_SUPPORT)
  return()
endif()
message("Building roshambo2 wrapper.")

if(NOT DEFINED ROSHAMBOSHAPE_DIR)
  set(ROSHAMBOSHAPE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/roshambo2")
  set(fileToCheck "${ROSHAMBOSHAPE_DIR}/roshambo2/backends/cpp_src/cpp_helper_functions.cpp")
  message("fileToCheck : ${fileToCheck}")
  set(needDownload "TRUE")
  if(EXISTS "${fileToCheck}")
    set(needDownload "FALSE")
  endif()
else()
  set(needDownload "FALSE")
endif()

message("Need download of roshambo2 : ${needDownload}")

if(needDownload)
  set(ROSHAMBOSHAPE_BASE "DecoupleCPPInterface.tar.gz")
  if(NOT DEFINED ROSHAMBOSHAPE_URL)
    set(ROSHAMBOSHAPE_URL "https://github.com/davidacosgrove/roshambo2/archive/${ROSHAMBOSHAPE_BASE}")
  endif()
  if(NOT DEFINED ROSHAMBOSHAPE_MD5SUM)
    set(ROSHAMBOSHAPE_MD5SUM "6b96a114704298cad266ff359045a73d")
  endif()
  downloadAndCheckMD5(${ROSHAMBOSHAPE_URL} "${CMAKE_CURRENT_SOURCE_DIR}/${ROSHAMBOSHAPE_BASE}" ${ROSHAMBOSHAPE_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf
    ${CMAKE_CURRENT_SOURCE_DIR}/${ROSHAMBOSHAPE_BASE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  message("${tar_dirname} and ${ROSHAMBOSHAPE_DIR}")
  file(GLOB tar_dirname ${CMAKE_CURRENT_SOURCE_DIR}/roshambo2-DecoupleCPPInterface*)
  execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${tar_dirname}
       ${ROSHAMBOSHAPE_DIR}
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

rdkit_library(Roshambo2Shape Roshambo2Shape.cpp ShapeInput.cpp ShapeOverlayOptions.cpp
        ./roshambo2/roshambo2/backends/cpp_src/cpp_helper_functions.cpp SHARED
        LINK_LIBRARIES SmilesParse SubstructMatch MolTransforms)
target_compile_definitions(Roshambo2Shape PRIVATE RDKIT_PUBCHEMSHAPE_BUILD)

rdkit_headers(Roshambo2Shape.hpp ShapeInput.h ShapeOverlayOptions.h DEST GraphMol)

rdkit_catch_test(roshambo2_test catch_tests.cpp LINK_LIBRARIES Roshambo2Shape FileParsers
        MolAlign MolTransforms)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()
