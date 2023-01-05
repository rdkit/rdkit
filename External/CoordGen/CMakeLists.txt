add_custom_target(coordgen_support ALL)
include(RDKitUtils)

if(RDK_BUILD_MAEPARSER_SUPPORT)
  if((NOT MSVC) OR RDK_INSTALL_DLLS_MSVC)
    add_definitions(-DIN_MAEPARSER)
  endif()

  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set (link_iostreams ${Boost_LIBRARIES})
  if (NOT Boost_USE_STATIC_LIBS)
    add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
  endif()

  if (WIN32)
    set( T_LIBS ${Boost_LIBRARIES})
    find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS zlib)
    set( Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})

  else()
    if(Boost_USE_STATIC_LIBS)
      # when we're doing static linkage of boost
      # to also link against zlib (due to iostreams)
      find_package(ZLIB)
    endif()
  endif()

  find_package(maeparser MODULE)

  if(MAEPARSER_FORCE_BUILD OR (NOT maeparser_FOUND))
    if(NOT DEFINED MAEPARSER_DIR)
      set(MAEPARSER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser")
    endif()

    if(NOT EXISTS "${MAEPARSER_DIR}/MaeParser.hpp")
        set(RELEASE_NO "1.2.4")
        set(MD5 "42a0765f9ddac81d588f1dc5afc87ed4")
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

    rdkit_library(maeparser ${MAESOURCES} LINK_LIBRARIES ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} SHARED)
    install(TARGETS maeparser DESTINATION ${RDKit_LibDir})
    set(maeparser_LIBRARIES maeparser ${ZLIB_LIBRARIES})


  endif(MAEPARSER_FORCE_BUILD OR (NOT maeparser_FOUND))

  include_directories(${maeparser_INCLUDE_DIRS})

  set(RDK_MAEPARSER_LIBS ${maeparser_LIBRARIES} ${Boost_LIBRARIES}
      CACHE STRING "the external libraries" FORCE)

else (RDK_BUILD_MAEPARSER_SUPPORT)

  set(RDK_MAEPARSER_LIBS CACHE STRING "the external libraries" FORCE)

endif(RDK_BUILD_MAEPARSER_SUPPORT)

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
        set(RELEASE_NO "3.0.1")
        set(MD5 "60951499eee3f55778ce8fe3ae190fbc")
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
    rdkit_library(coordgen ${CGSOURCES} SHARED )
    install(TARGETS coordgen DESTINATION ${RDKit_LibDir})
    set(coordgen_LIBRARIES coordgen)

  endif(COORDGEN_FORCE_BUILD OR (NOT coordgen_FOUND))

  include_directories(${coordgen_INCLUDE_DIRS})

  set(RDK_COORDGEN_LIBS MolAlign ${coordgen_LIBRARIES} ${RDK_MAEPARSER_LIBS} ${Boost_LIBRARIES}
      CACHE STRING "the external libraries" FORCE)
  rdkit_headers(CoordGen.h DEST GraphMol)

  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  rdkit_test(testCoordGen test.cpp
    LINK_LIBRARIES
    ${RDK_COORDGEN_LIBS} Depictor ChemTransforms
    FileParsers SmilesParse SubstructMatch GraphMol
    RDGeneral DataStructs RDGeneral RDGeometryLib
    ${RDKit_THREAD_LIBS})

else (RDK_BUILD_COORDGEN_SUPPORT)

  set(RDK_COORDGEN_LIBS CACHE STRING "the external libraries" FORCE)

endif(RDK_BUILD_COORDGEN_SUPPORT)
