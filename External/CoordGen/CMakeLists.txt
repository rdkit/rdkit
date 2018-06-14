add_custom_target(coordgen_support ALL)

if(RDK_BUILD_COORDGEN_SUPPORT)
    add_definitions(-DIN_MAEPARSER)
    add_definitions(-DIN_COORDGEN)
    include(RDKitUtils)
    if(NOT DEFINED MAEPARSER_DIR)
      set(MAEPARSER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser")
    endif()
    if(NOT EXISTS "${MAEPARSER_DIR}/MaeParser.hpp")
        set(RELEASE_NO "1.0.0")
        downloadAndCheckMD5("https://github.com/schrodinger/maeparser/archive/${RELEASE_NO}.tar.gz"
              "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
              "770d941264dde4825aaa505e59642e51")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "maeparser-${RELEASE_NO}" "${MAEPARSER_DIR}")
        patchCoordGenMaeExportHeaders("MAEPARSER" "${MAEPARSER_DIR}/MaeParserConfig.hpp")
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
        set(RELEASE_NO "1.1")
        set(SHA "ede3191bf1c1dde53f70ca086014ff9ec7768c55")
        downloadAndCheckMD5("https://github.com/schrodinger/coordgenlibs/archive/v${RELEASE_NO}.tar.gz"
              "${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz"
              "493cf2eab7d4ab6242e25da112831798")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/master.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "coordgenlibs-${RELEASE_NO}" "${COORDGEN_DIR}")
        patchCoordGenMaeExportHeaders("COORDGEN" "${COORDGEN_DIR}/CoordgenConfig.hpp")
    else()
      message("-- Found coordgenlibs source in ${COORDGEN_DIR}")
    endif()

    file(GLOB CGSOURCES "${COORDGEN_DIR}/*.cpp")
    rdkit_library(coordgenlib ${CGSOURCES} SHARED LINK_LIBRARIES maeparser)
    install(TARGETS coordgenlib DESTINATION ${RDKit_LibDir})

    set(MAE_FILES ${COORDGEN_DIR}/templates.mae )
    install(FILES ${MAE_FILES}
            DESTINATION ${RDKit_ShareDir}/Data
            COMPONENT data )

    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    set(RDK_COORDGEN_LIBS MolAlign coordgenlib maeparser CACHE STRING "the external libraries" FORCE)
    rdkit_headers(CoordGen.h DEST GraphMol)

    if(RDK_BUILD_PYTHON_WRAPPERS)
      add_subdirectory(Wrap)
    endif(RDK_BUILD_PYTHON_WRAPPERS)

    rdkit_test(testCoordGen test.cpp
      LINK_LIBRARIES
      ${RDK_COORDGEN_LIBS} Depictor
      FileParsers SmilesParse SubstructMatch GraphMol
      RDGeneral DataStructs RDGeneral RDGeometryLib
      ${RDKit_THREAD_LIBS})
else (RDK_BUILD_COORDGEN_SUPPORT)
    set(RDK_COORDGEN_LIBS  CACHE STRING "the external libraries" FORCE)
endif(RDK_BUILD_COORDGEN_SUPPORT)
