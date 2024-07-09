add_custom_target(coordgen_support ALL)
include(RDKitUtils)

if(RDK_BUILD_MAEPARSER_SUPPORT OR RDK_BUILD_COORDGEN_SUPPORT)
  if((NOT MSVC) OR RDK_INSTALL_DLLS_MSVC)
    add_definitions(-DIN_MAEPARSER)
  endif()

  find_package(maeparser MODULE)

  if(MAEPARSER_FORCE_BUILD OR (NOT maeparser_FOUND))
    if(NOT DEFINED MAEPARSER_DIR)
      set(MAEPARSER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser")
    endif()

    if(NOT EXISTS "${MAEPARSER_DIR}/MaeParser.hpp")
        set(RELEASE_NO "1.3.1")
        set(MD5 "cfa40e29366f4b413e4ec15f959ee139")
        downloadAndCheckMD5("https://github.com/schrodinger/maeparser/archive/v${RELEASE_NO}.tar.gz"
              "${CMAKE_CURRENT_SOURCE_DIR}/maeparser-v${RELEASE_NO}.tar.gz" ${MD5})
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/maeparser-v${RELEASE_NO}.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "maeparser-${RELEASE_NO}" "${MAEPARSER_DIR}")
        patchCoordGenMaeExportHeaders("MAEPARSER" "${MAEPARSER_DIR}/MaeParserConfig.hpp")

    else()
      message("-- Found MAEParser source in ${MAEPARSER_DIR}")
    endif()

    set(maeparser_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        CACHE STRING "MaeParser Include Dir" FORCE)
    file(GLOB MAESOURCES "${MAEPARSER_DIR}/*.cpp")

     if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-copy")
     endif()

    rdkit_library(maeparser ${MAESOURCES} SHARED)
    install(TARGETS maeparser DESTINATION ${RDKit_LibDir})
    set(maeparser_LIBRARIES maeparser)


  endif(MAEPARSER_FORCE_BUILD OR (NOT maeparser_FOUND))

  include_directories(${maeparser_INCLUDE_DIRS})

  set(RDK_MAEPARSER_LIBS ${maeparser_LIBRARIES}
      CACHE STRING "the external libraries" FORCE)

else ()

  set(RDK_MAEPARSER_LIBS CACHE STRING "the external libraries" FORCE)

endif()

if(RDK_BUILD_COORDGEN_SUPPORT)
  if(MSVC AND (NOT RDK_INSTALL_DLLS_MSVC))
    add_definitions(-DSTATIC_COORDGEN)
  else(MSVC AND (NOT RDK_INSTALL_DLLS_MSVC))
    add_definitions(-DIN_COORDGEN)
  endif(MSVC AND (NOT RDK_INSTALL_DLLS_MSVC))

  find_package(coordgen MODULE)
  if(COORDGEN_FORCE_BUILD OR (NOT coordgen_FOUND))
    if(NOT DEFINED COORDGEN_DIR)
      set(COORDGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coordgen")
    endif()

    if(NOT EXISTS "${COORDGEN_DIR}/sketcherMinimizer.h")
        set(RELEASE_NO "3.0.2")
        set(MD5 "bc9dabbbb2b99ae426f3a1bf16c4d3cc")
        downloadAndCheckMD5("https://github.com/schrodinger/coordgenlibs/archive/v${RELEASE_NO}.tar.gz"
              "${CMAKE_CURRENT_SOURCE_DIR}/coordgenlibs-${RELEASE_NO}.tar.gz" ${MD5})
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/coordgenlibs-${RELEASE_NO}.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        file(RENAME "coordgenlibs-${RELEASE_NO}" "${COORDGEN_DIR}")
    else()
      message("-- Found coordgenlibs source in ${COORDGEN_DIR}")
    endif()

    set(coordgen_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        CACHE STRING "CoordGen Include File" FORCE)
    file(GLOB CGSOURCES "${COORDGEN_DIR}/*.cpp")
    rdkit_library(coordgen ${CGSOURCES} SHARED)
if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(coordgen PUBLIC -Wno-unused-but-set-variable)
if(RDK_INSTALL_STATIC_LIBS)
    target_compile_options(coordgen_static PUBLIC -Wno-unused-but-set-variable)
endif()
endif()
    install(TARGETS coordgen DESTINATION ${RDKit_LibDir})
    set(RDK_COORDGEN_LIBS coordgen CACHE STRING "the external libraries" FORCE)

  elseif(coordgen_FOUND)

    set(RDK_COORDGEN_LIBS "${coordgen_LIBRARIES}" CACHE STRING "the external libraries" FORCE)

  endif(COORDGEN_FORCE_BUILD OR(NOT coordgen_FOUND))

  include_directories(${coordgen_INCLUDE_DIRS})

  rdkit_headers(CoordGen.h DEST GraphMol)

  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  rdkit_test(testCoordGen test.cpp
    LINK_LIBRARIES
    "${RDK_COORDGEN_LIBS}" Depictor ChemTransforms
    FileParsers SmilesParse SubstructMatch GraphMol
    RDGeneral DataStructs RDGeneral RDGeometryLib
    ${RDKit_THREAD_LIBS})

else (RDK_BUILD_COORDGEN_SUPPORT)

  set(RDK_COORDGEN_LIBS CACHE STRING "the external libraries" FORCE)

endif(RDK_BUILD_COORDGEN_SUPPORT)
