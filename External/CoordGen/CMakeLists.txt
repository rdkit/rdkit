add_custom_target(coordgen_support ALL)

if(RDK_BUILD_COORDGEN_SUPPORT)
    if(NOT DEFINED MAEPARSER_DIR)
      set(MAEPARSER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser")
    endif()
    if(NOT EXISTS "${MAEPARSER_DIR}/MaeParser.hpp")
        set(SHA "83368293dcc0eb07562dadfb7728b8d18d23a6cb")
        downloadAndCheckMD5("https://codeload.github.com/schrodinger/maeparser/tar.gz/${SHA}"
              "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
              "32c0c3b315bba49fbf4c41a07aa58528")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "maeparser-${SHA}" "${MAEPARSER_DIR}")
    else()
      message("-- Found MAEParser source in ${MAEPARSER_DIR}")
    endif()

    file(GLOB MAESOURCES "${MAEPARSER_DIR}/*.cpp")
    rdkit_library(maeparser ${MAESOURCES} SHARED )
    install(TARGETS maeparser DESTINATION ${RDKit_LibDir})

    if(NOT DEFINED COORDGEN_DIR)
      set(COORDGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coordgenlibs")
    endif()
    if(NOT EXISTS "${COORDGEN_DIR}/sketcherMinimizer.h")
        set(SHA "68fb2083d91a3f1de734e2c4f2f7f2d98e21634d")
        downloadAndCheckMD5("https://codeload.github.com/schrodinger/coordgenlibs/tar.gz/${SHA}"
              "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
              "3e8b644c9a52333cef506724db51ffc1")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "coordgenlibs-${SHA}" "${COORDGEN_DIR}")
    else()
      message("-- Found coordgenlibs source in ${COORDGEN_DIR}")
    endif()

    file(GLOB CGSOURCES "${COORDGEN_DIR}/*.cpp")
    rdkit_library(coordgenlib ${CGSOURCES} SHARED LINK_LIBRARIES maeparser)
    install(TARGETS coordgenlib DESTINATION ${RDKit_LibDir})

    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    set(COORDGEN_LIBS coordgenlib maeparser )
  #rdkit_library(RDCoordGenLib coordgen.cpp SHARED LINK_LIBRARIES ${INCHI_LIBRARIES}
  #              GraphMol RDGeneral Depictor SubstructMatch SmilesParse
  #              ${RDKit_THREAD_LIBS})
  rdkit_headers(CoordGen.h DEST GraphMol)
  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  rdkit_test(testCoordGen test.cpp
    LINK_LIBRARIES
    ${COORDGEN_LIBS}
    FileParsers SmilesParse SubstructMatch GraphMol
    RDGeneral DataStructs RDGeneral RDGeometryLib
    ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_COORDGEN_SUPPORT)
