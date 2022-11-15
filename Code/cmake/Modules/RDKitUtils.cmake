include(BoostUtils)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac OS X specific code
  set(RDKit_VERSION "${RDKit_Year}.${RDKit_Month}")
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RDKit_VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(RDKit_VERSION "${RDKit_VERSION}.${RDKit_Revision}${RDKit_RevisionModifier}")
set(RDKit_RELEASENAME "${RDKit_Year}.${RDKit_Month}.${RDKit_Revision}${RDKit_RevisionModifier}")


set(compilerID "${CMAKE_CXX_COMPILER_ID}")
set(systemAttribute "")
if(MINGW)
  set(systemAttribute "MINGW")
endif(MINGW)
if(UNIX)
  set(systemAttribute "UNIX")
endif(UNIX)
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(bit3264 "32-bit")
else()
  set(bit3264 "64-bit")
endif()
set(RDKit_BUILDNAME "${CMAKE_SYSTEM_NAME}|${CMAKE_SYSTEM_VERSION}|${systemAttribute}|${compilerID}|${bit3264}")
set(RDKit_EXPORTED_TARGETS rdkit-targets)


macro(rdkit_library)
  PARSE_ARGUMENTS(RDKLIB
    "LINK_LIBRARIES;DEST"
    "SHARED"
    ${ARGN})
  CAR(RDKLIB_NAME ${RDKLIB_DEFAULT_ARGS})
  CDR(RDKLIB_SOURCES ${RDKLIB_DEFAULT_ARGS})
  if((MSVC AND (NOT RDK_INSTALL_DLLS_MSVC)) OR (WIN32 AND RDK_INSTALL_STATIC_LIBS))
    add_library(${RDKLIB_NAME} ${RDKLIB_SOURCES})
    target_link_libraries(${RDKLIB_NAME} PUBLIC rdkit_base)
    if(RDK_INSTALL_DEV_COMPONENT)
      INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
              DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
              COMPONENT dev )
    endif(RDK_INSTALL_DEV_COMPONENT)
  else()
    # we're going to always build in shared mode since we
    # need exceptions to be (correctly) catchable across
    # boundaries. As of now (June 2010), this doesn't work
    # with g++ unless libraries are shared.
    add_library(${RDKLIB_NAME} SHARED ${RDKLIB_SOURCES})
    target_link_libraries(${RDKLIB_NAME} PUBLIC rdkit_base)
    INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
            DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
            COMPONENT runtime )
    if(RDK_INSTALL_STATIC_LIBS)
      add_library(${RDKLIB_NAME}_static ${RDKLIB_SOURCES})

      set(skipNext FALSE)
      foreach(linkLib ${RDKLIB_LINK_LIBRARIES})
        if(skipNext)
          set(skipNext FALSE)
          continue()
        endif()
        if(TARGET "${linkLib}")
          get_target_property(linkLib_IMPORTED "${linkLib}" IMPORTED)
          if (linkLib_IMPORTED)
            # linkLib is an imported target: use it as-is
            target_link_libraries(${RDKLIB_NAME}_static PUBLIC "${linkLib}")
            continue()
          endif()
        elseif(EXISTS "${linkLib}")
          # linkLib is a file, so keep it as-is
          target_link_libraries(${RDKLIB_NAME}_static PUBLIC "${linkLib}")
          continue()
        # cmake prepends the special keywords debug, optimized, general
        # before the library name depending on whether they should be
        # linked in Debug, Release or generic builds. Therefore we need
        # to skip those, and also skip the library that follows if it
        # is not relevant for the current build type
        elseif ("${linkLib}" STREQUAL "debug")
          if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(skipNext TRUE)
          endif()
          continue()
        elseif ("${linkLib}" STREQUAL "optimized")
          if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(skipNext TRUE)
          endif()
          continue()
        elseif ("${linkLib}" STREQUAL "general")
          continue()
        endif()

        # We haven't seen linkLib yet. This probably means it is a target
        # we will be creating at some point (if not, then we are missing a find_package).
        # Add the "_static" suffix to link against its static variant
        target_link_libraries(${RDKLIB_NAME}_static PUBLIC "${linkLib}_static")
      endforeach()
      target_link_libraries(${RDKLIB_NAME}_static PUBLIC rdkit_base)
      if(RDK_INSTALL_DEV_COMPONENT)
        INSTALL(TARGETS ${RDKLIB_NAME}_static EXPORT ${RDKit_EXPORTED_TARGETS}
                DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
                COMPONENT dev )
      endif(RDK_INSTALL_DEV_COMPONENT)
      set_target_properties(${RDKLIB_NAME}_static PROPERTIES
                            OUTPUT_NAME "RDKit${RDKLIB_NAME}_static")

    endif(RDK_INSTALL_STATIC_LIBS)
  endif()
  IF(RDKLIB_LINK_LIBRARIES)
    target_link_libraries(${RDKLIB_NAME} PUBLIC ${RDKLIB_LINK_LIBRARIES})
  ENDIF(RDKLIB_LINK_LIBRARIES)
  if((NOT MSVC) OR RDK_INSTALL_DLLS_MSVC)
    set_target_properties(${RDKLIB_NAME} PROPERTIES
                          OUTPUT_NAME "RDKit${RDKLIB_NAME}"
                          VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}.${RDKit_Revision}"
                          VERSION ${RDKit_VERSION}
                          SOVERSION ${RDKit_ABI} )
  endif()
  set_target_properties(${RDKLIB_NAME} PROPERTIES
                        ARCHIVE_OUTPUT_DIRECTORY ${RDK_ARCHIVE_OUTPUT_DIRECTORY}
                        RUNTIME_OUTPUT_DIRECTORY ${RDK_RUNTIME_OUTPUT_DIRECTORY}
                        LIBRARY_OUTPUT_DIRECTORY ${RDK_LIBRARY_OUTPUT_DIRECTORY})
endmacro(rdkit_library)

macro(rdkit_headers)
  if (NOT RDK_INSTALL_INTREE)
    PARSE_ARGUMENTS(RDKHDR
      "DEST"
      ""
      ${ARGN})
    # RDKHDR_DEFAULT_ARGS -> RDKHDR_DEST
    if(RDK_INSTALL_DEV_COMPONENT)
      install(FILES ${RDKHDR_DEFAULT_ARGS}
              DESTINATION ${RDKit_HdrDir}/${RDKHDR_DEST}
              COMPONENT dev )
    endif(RDK_INSTALL_DEV_COMPONENT)
  endif(NOT RDK_INSTALL_INTREE)
endmacro(rdkit_headers)

macro(rdkit_python_extension)
  PARSE_ARGUMENTS(RDKPY
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(RDKPY_NAME ${RDKPY_DEFAULT_ARGS})
  CDR(RDKPY_SOURCES ${RDKPY_DEFAULT_ARGS})
  if(RDK_BUILD_PYTHON_WRAPPERS)
    PYTHON_ADD_MODULE(${RDKPY_NAME} ${RDKPY_SOURCES})
    set_target_properties(${RDKPY_NAME} PROPERTIES PREFIX "")

    if(WIN32)
      set_target_properties(${RDKPY_NAME} PROPERTIES SUFFIX ".pyd"
                           LIBRARY_OUTPUT_DIRECTORY
                           ${RDK_PYTHON_OUTPUT_DIRECTORY}/${RDKPY_DEST})
    else(WIN32)
        set_target_properties(${RDKPY_NAME} PROPERTIES
                              LIBRARY_OUTPUT_DIRECTORY
                              ${RDK_PYTHON_OUTPUT_DIRECTORY}/${RDKPY_DEST})
    endif(WIN32)

    target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES}
                          RDBoost rdkit_py_base rdkit_base )
    if("${PYTHON_LDSHARED}" STREQUAL "")
    else()
      set_target_properties(${RDKPY_NAME} PROPERTIES LINK_FLAGS ${PYTHON_LDSHARED})
    endif()

    INSTALL(TARGETS ${RDKPY_NAME}
            LIBRARY DESTINATION ${RDKit_PythonDir}/${RDKPY_DEST} COMPONENT python)
  endif(RDK_BUILD_PYTHON_WRAPPERS)
endmacro(rdkit_python_extension)

macro(rdkit_test)
  PARSE_ARGUMENTS(RDKTEST
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(RDKTEST_NAME ${RDKTEST_DEFAULT_ARGS})
  CDR(RDKTEST_SOURCES ${RDKTEST_DEFAULT_ARGS})
  if(RDK_BUILD_CPP_TESTS)
    add_executable(${RDKTEST_NAME} ${RDKTEST_SOURCES})
    target_link_libraries(${RDKTEST_NAME} ${RDKTEST_LINK_LIBRARIES})
    add_test(${RDKTEST_NAME} ${EXECUTABLE_OUTPUT_PATH}/${RDKTEST_NAME})
  endif(RDK_BUILD_CPP_TESTS)
endmacro(rdkit_test)

macro(rdkit_catch_test)
  PARSE_ARGUMENTS(RDKTEST
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(RDKTEST_NAME ${RDKTEST_DEFAULT_ARGS})
  CDR(RDKTEST_SOURCES ${RDKTEST_DEFAULT_ARGS})
  if(RDK_BUILD_CPP_TESTS)
    add_executable(${RDKTEST_NAME} ${RDKTEST_SOURCES})
    target_include_directories(${RDKTEST_NAME} PRIVATE ${CATCH_INCLUDE_DIR})
    target_link_libraries(${RDKTEST_NAME} rdkitCatch ${RDKTEST_LINK_LIBRARIES})
    add_test(${RDKTEST_NAME} ${EXECUTABLE_OUTPUT_PATH}/${RDKTEST_NAME})
    #ParseAndAddCatchTests(${RDKTEST_NAME})
    add_dependencies(${RDKTEST_NAME} catch)
  endif(RDK_BUILD_CPP_TESTS)
endmacro(rdkit_catch_test)

macro(add_pytest)
  PARSE_ARGUMENTS(PYTEST
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(PYTEST_NAME ${PYTEST_DEFAULT_ARGS})
  CDR(PYTEST_SOURCES ${PYTEST_DEFAULT_ARGS})
  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_test(${PYTEST_NAME}  ${PYTHON_EXECUTABLE}
             ${PYTEST_SOURCES})
    SET(RDKIT_PYTEST_CACHE "${PYTEST_NAME};${RDKIT_PYTEST_CACHE}" CACHE INTERNAL "Global list of python tests")
  endif(RDK_BUILD_PYTHON_WRAPPERS)
endmacro(add_pytest)

function(add_jupytertest testname workingdir notebook)
  if(RDK_BUILD_PYTHON_WRAPPERS AND RDK_NBVAL_AVAILABLE)
    add_test(NAME ${testname}  COMMAND ${PYTHON_EXECUTABLE} -m pytest --nbval ${notebook}
       WORKING_DIRECTORY ${workingdir} )
    SET(RDKIT_JUPYTERTEST_CACHE "${testname};${RDKIT_JUPYTERTEST_CACHE}" CACHE INTERNAL "Global list of jupyter tests")
  endif()
endfunction(add_jupytertest)

function(computeMD5 target md5chksum)
  execute_process(COMMAND ${CMAKE_COMMAND} -E md5sum ${target} OUTPUT_VARIABLE md5list)
  string(REGEX REPLACE "([a-z0-9]+)" "\\1;" md5list "${md5list}")
  list(GET md5list 0 md5)
  set(${md5chksum} ${md5} PARENT_SCOPE)
endfunction(computeMD5)

function(downloadAndCheckMD5 url target md5chksum)
  if ("${url}" STREQUAL "")
    return()
  endif()
  get_filename_component(targetDir ${target} PATH)
  message("Downloading ${url}...")
  file(DOWNLOAD "${url}" "${target}"
    STATUS status)
  # CMake < 2.8.10 does not seem to support HTTPS out of the box
  # and since SourceForge redirects to HTTPS, the CMake download fails
  # so we try to use Powershell (Windows) or system curl (Unix, OS X) if available
  if (NOT status EQUAL 0)
    if(WIN32)
      execute_process(COMMAND powershell -Command "(New-Object Net.WebClient).DownloadFile('${url}', '${target}')")
    else(WIN32)
      execute_process(COMMAND curl -L "${url}" -o ${target} WORKING_DIRECTORY ${targetDir})
    endif(WIN32)
  endif()
  if (NOT EXISTS ${target})
    MESSAGE(FATAL_ERROR "The download of ${url} failed.")
  endif()
  if (NOT "${md5chksum}" STREQUAL "")
    computeMD5(${target} md5)
    if (NOT md5 STREQUAL ${md5chksum})
      MESSAGE(FATAL_ERROR "The md5 checksum for ${target} is incorrect; expected: ${md5chksum}, found: ${md5}")
    endif()
  endif()
endfunction(downloadAndCheckMD5)

function(createExportTestHeaders)
  file(GLOB_RECURSE cmakeLists LIST_DIRECTORIES false
       ${CMAKE_SOURCE_DIR}/CMakeLists.txt)
  set(exportLibs "")
  foreach(cmakeList ${cmakeLists})
    file(STRINGS ${cmakeList} rdkitLibraryItems REGEX "rdkit_library[ ]*\\([ ]*[^ ]+.*$")
    if (NOT "${rdkitLibraryItems}" STREQUAL "")
      foreach (rdkitLibrary ${rdkitLibraryItems})
        string(REGEX REPLACE "^[ ]*rdkit_library[ ]*\\([ ]*([^ ]+).*$" "\\1" libName "${rdkitLibrary}")
        list(APPEND exportLibs "${libName}")
      endforeach()
    endif()
  endforeach()
  list(REMOVE_DUPLICATES exportLibs)
  list(SORT exportLibs)
  set(exportPath "Code/RDGeneral/export.h")
  file(WRITE "${CMAKE_BINARY_DIR}/${exportPath}.tmp"
    "// auto-generated export definition header\n"
    "#pragma once\n"
    "#include <RDGeneral/RDExportMacros.h>\n")
  set(testPath "Code/RDGeneral/test.h")
  file(WRITE "${CMAKE_BINARY_DIR}/${testPath}.tmp"
    "// auto-generated header to be imported in all cpp tests\n"
    "#pragma once\n")
  foreach(exportLib ${exportLibs})
    string(TOUPPER "${exportLib}" exportLib)
    file(APPEND "${CMAKE_BINARY_DIR}/${exportPath}.tmp"
      "\n"
      "// RDKIT_${exportLib}_EXPORT definitions\n"
      "#ifdef RDKIT_${exportLib}_BUILD\n"
      "#define RDKIT_${exportLib}_EXPORT RDKIT_EXPORT_API\n"
      "#else\n"
      "#define RDKIT_${exportLib}_EXPORT RDKIT_IMPORT_API\n"
      "#endif\n"
      "// RDKIT_${exportLib}_EXPORT end definitions\n")
    file(APPEND "${CMAKE_BINARY_DIR}/${testPath}.tmp"
      "\n"
      "#ifdef RDKIT_${exportLib}_BUILD\n"
      "#undef RDKIT_${exportLib}_BUILD\n"
      "#endif\n")
  endforeach()
  file(APPEND "${CMAKE_BINARY_DIR}/${exportPath}.tmp"
  "\n"
  "/*\n"
  " * Do not dll export/import to export queries (it will mess up with the\n"
  " * templates), but make sure it is visible for *nix\n"
  " */\n"
  "// RDKIT_QUERY_EXPORT definitions\n"
  "#if defined(RDKIT_DYN_LINK) && defined(WIN32) && defined(BOOST_HAS_DECLSPEC)\n"
  "#define RDKIT_QUERY_EXPORT\n"
  "#else\n"
  "#define RDKIT_QUERY_EXPORT RDKIT_GRAPHMOL_EXPORT\n"
  "#endif\n"
  "// RDKIT_QUERY_EXPORT end definitions\n")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_BINARY_DIR}/${exportPath}.tmp" "${CMAKE_SOURCE_DIR}/${exportPath}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_BINARY_DIR}/${testPath}.tmp" "${CMAKE_SOURCE_DIR}/${testPath}")
endfunction(createExportTestHeaders)

function(patchCoordGenMaeExportHeaders keyword path)
  file(APPEND "${path}"
    "// appended by CMake patchCoordGenMaeExportHeaders\n"
    "#if !defined(RDKIT_DYN_LINK) || defined(SWIG)\n"
    "#ifdef EXPORT_${keyword}\n"
    "#undef EXPORT_${keyword}\n"
    "#endif\n"
    "#define EXPORT_${keyword}\n"
    "#endif\n"
    "#ifndef SWIG\n"
    "#ifdef _MSC_VER\n"
    "#pragma warning(disable:4251)\n"
    "#pragma warning(disable:4275)\n"
    "#endif\n"
    "#endif\n")
endfunction(patchCoordGenMaeExportHeaders)
