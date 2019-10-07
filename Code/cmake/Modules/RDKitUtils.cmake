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
if(CMAKE_SIZEOF_VOID_P MATCHES 4)
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
  if(MSVC AND (NOT RDK_INSTALL_DLLS_MSVC))
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

      foreach(linkLib ${RDKLIB_LINK_LIBRARIES})
        if(${linkLib} MATCHES "^(Boost)|(Thread)|(boost)|^(optimized)|^(debug)|(libz)")
          set(rdk_static_link_libraries "${rdk_static_link_libraries}${linkLib};")
        else()
          set(rdk_static_link_libraries "${rdk_static_link_libraries}${linkLib}_static;")
        endif()
      endforeach()
      target_link_libraries(${RDKLIB_NAME}_static PUBLIC ${rdk_static_link_libraries})
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

    if(WIN32 OR "${Py_ENABLE_SHARED}" STREQUAL "1")
      target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES}
                            ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} )
    else()
      target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES}
                            ${Boost_LIBRARIES} )
      if("${PYTHON_LDSHARED}" STREQUAL "")
      else()
        set_target_properties(${RDKPY_NAME} PROPERTIES LINK_FLAGS ${PYTHON_LDSHARED})
      endif()
    endif()

    if(RDK_INSTALL_INTREE)
      INSTALL(TARGETS ${RDKPY_NAME}
              LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/rdkit/${RDKPY_DEST} COMPONENT python)
    else(RDK_INSTALL_INTREE)
      file(MAKE_DIRECTORY ${RDK_PYTHON_BINARY_DIR}/rdkit/${RDKPY_DEST})

      # Copy built library into Python package tree.
      add_custom_command(TARGET ${RDKPY_NAME} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${RDKPY_NAME}> ${RDK_PYTHON_BINARY_DIR}/rdkit/${RDKPY_DEST}
      )
      add_dependencies(python_dist ${RDKPY_NAME})
      # We will build .whl package with these libs at __build time__, which means they already must have
      # correct RPATHs.
      set_target_properties(${RDKPY_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
      set_target_properties(${RDKPY_NAME} PROPERTIES BUILD_WITH_INSTALL_NAME_DIR ON)
    endif(RDK_INSTALL_INTREE)
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
    target_link_libraries(${RDKTEST_NAME} ${RDKTEST_LINK_LIBRARIES})
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

function(computeMD5 target md5chksum)
  execute_process(COMMAND ${CMAKE_COMMAND} -E md5sum ${target} OUTPUT_VARIABLE md5list)
  string(REGEX REPLACE "([a-z0-9]+)" "\\1;" md5list "${md5list}")
  list(GET md5list 0 md5)
  set(${md5chksum} ${md5} PARENT_SCOPE)
endfunction(computeMD5)

function(downloadAndCheckMD5 url target md5chksum)
  if (NOT ${url} EQUAL "")
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
    if (NOT ${md5chksum} EQUAL "")
      computeMD5(${target} md5)
      if (NOT md5 STREQUAL ${md5chksum})
        MESSAGE(FATAL_ERROR "The md5 checksum for ${target} is incorrect; expected: ${md5chksum}, found: ${md5}")
      endif()
    endif()
  endif()
endfunction(downloadAndCheckMD5)

function(overwriteIfChanged src dest)
  set(overwrite TRUE)
  if (EXISTS "${dest}")
    computeMD5("${dest}" destMD5)
    computeMD5("${src}" srcMD5)
    if (${destMD5} STREQUAL ${srcMD5})
      set(overwrite FALSE)
    endif()
  endif()
  if (overwrite)
    get_filename_component(destDir "${dest}" DIRECTORY)
    file(COPY "${src}" DESTINATION "${destDir}")
  endif()
endfunction(overwriteIfChanged)

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
  file(WRITE "${CMAKE_BINARY_DIR}/${exportPath}"
    "// auto-generated __declspec definition header\n"
    "#pragma once\n"
    "#ifndef SWIG\n"
    "#ifdef _MSC_VER\n"
    "#pragma warning(disable:4251)\n"
    "#pragma warning(disable:4275)\n"
    "#endif\n"
    "\n"
    "#include <boost/config.hpp>\n"
    "#endif\n")
  set(testPath "Code/RDGeneral/test.h")
  file(WRITE "${CMAKE_BINARY_DIR}/${testPath}"
    "// auto-generated header to be imported in all cpp tests\n"
    "#pragma once\n")
  foreach(exportLib ${exportLibs})
    string(TOUPPER "${exportLib}" exportLib)
    file(APPEND "${CMAKE_BINARY_DIR}/${exportPath}"
      "\n"
      "// RDKIT_${exportLib}_EXPORT definitions\n"
      "#if defined(BOOST_HAS_DECLSPEC) && defined(RDKIT_DYN_LINK) && !defined(SWIG)\n"
      "#ifdef RDKIT_${exportLib}_BUILD\n"
      "#define RDKIT_${exportLib}_EXPORT __declspec(dllexport)\n"
      "#else\n"
      "#define RDKIT_${exportLib}_EXPORT __declspec(dllimport)\n"
      "#endif\n"
      "#endif\n"
      "#ifndef RDKIT_${exportLib}_EXPORT\n"
      "#define RDKIT_${exportLib}_EXPORT\n"
      "#endif\n"
      "// RDKIT_${exportLib}_EXPORT end definitions\n")
    file(APPEND "${CMAKE_BINARY_DIR}/${testPath}"
      "\n"
      "#ifdef RDKIT_${exportLib}_BUILD\n"
      "#undef RDKIT_${exportLib}_BUILD\n"
      "#endif\n")
  endforeach()
  overwriteIfChanged("${CMAKE_BINARY_DIR}/${exportPath}" "${CMAKE_SOURCE_DIR}/${exportPath}")
  overwriteIfChanged("${CMAKE_BINARY_DIR}/${testPath}" "${CMAKE_SOURCE_DIR}/${testPath}")
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
