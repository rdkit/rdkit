include(BoostUtils)

macro(rdkit_library)
  PARSE_ARGUMENTS(RDKLIB
    "LINK_LIBRARIES;DEST"
    "STATIC"
    ${ARGN})
  CAR(RDKLIB_NAME ${RDKLIB_DEFAULT_ARGS})
  CDR(RDKLIB_SOURCES ${RDKLIB_DEFAULT_ARGS})
  if(MSVC)
    add_library(${RDKLIB_NAME} ${RDKLIB_SOURCES})
  else(MSVC)
    IF(RDKLIB_STATIC)
      add_library(${RDKLIB_NAME} ${RDKLIB_SOURCES})
    ELSE(RDKLIB_STATIC)        
      add_library(${RDKLIB_NAME} SHARED ${RDKLIB_SOURCES})
      INSTALL(TARGETS ${RDKLIB_NAME} 
              DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST})
    ENDIF(RDKLIB_STATIC)        
    IF(RDKLIB_LINK_LIBRARIES)
      target_link_libraries(${RDKLIB_NAME} ${RDKLIB_LINK_LIBRARIES})
    ENDIF(RDKLIB_LINK_LIBRARIES)
  endif(MSVC)
endmacro(rdkit_library)
  
macro(rdkit_python_extension)
  PARSE_ARGUMENTS(RDKPY
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(RDKPY_NAME ${RDKPY_DEFAULT_ARGS})
  CDR(RDKPY_SOURCES ${RDKPY_DEFAULT_ARGS})
  PYTHON_ADD_MODULE(${RDKPY_NAME} ${RDKPY_SOURCES})
  set_target_properties(${RDKPY_NAME} PROPERTIES PREFIX "")
if(MSVC)
  set_target_properties(${RDKPY_NAME} PROPERTIES SUFFIX ".pyd")
endif(MSVC)  
  target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES} 
                        ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

  INSTALL(TARGETS ${RDKPY_NAME} 
          LIBRARY DESTINATION ${RDKit_PythonDir}/${RDKPY_DEST})
endmacro(rdkit_python_extension)

macro(rdkit_test)
  PARSE_ARGUMENTS(RDKTEST
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN})
  CAR(RDKTEST_NAME ${RDKTEST_DEFAULT_ARGS})
  CDR(RDKTEST_SOURCES ${RDKTEST_DEFAULT_ARGS})
  add_executable(${RDKTEST_NAME} ${RDKTEST_SOURCES})
  target_link_libraries(${RDKTEST_NAME} ${RDKTEST_LINK_LIBRARIES})
  add_test(${RDKTEST_NAME} ${EXECUTABLE_OUTPUT_PATH}/${RDKTEST_NAME})
endmacro(rdkit_test)

# ---------------
# downloaded from: http://www.cmake.org/Bug/file_download.php?file_id=2421&type=bug
# This additional function definition is needed to provide a workaround for
# CMake bug 9220.

# On debian testing (cmake 2.6.2), I get return code zero when calling 
# cmake the first time, but cmake crashes when running a second time
# as follows:
#
#  -- The Fortran compiler identification is unknown
#  CMake Error at /usr/share/cmake-2.6/Modules/CMakeFortranInformation.cmake:7 (GET_FILENAME_COMPONENT):
#    get_filename_component called with incorrect number of arguments
#  Call Stack (most recent call first):
#    CMakeLists.txt:3 (enable_language)
#
# My workaround is to invoke cmake twice.  If both return codes are zero, 
# it is safe to invoke ENABLE_LANGUAGE(Fortran OPTIONAL)

function(workaround_9220 language language_works)
  #message("DEBUG: language = ${language}")
  set(text
    "project(test NONE)
cmake_minimum_required(VERSION 2.6.0)
enable_language(${language} OPTIONAL)
"
    )
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/language_tests/${language})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language})
  file(WRITE ${CMAKE_BINARY_DIR}/language_tests/${language}/CMakeLists.txt
    ${text})
  execute_process(
    COMMAND ${CMAKE_COMMAND} .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language}
    RESULT_VARIABLE return_code
    OUTPUT_QUIET
    ERROR_QUIET
    )

  if(return_code EQUAL 0)
    # Second run
    execute_process (
      COMMAND ${CMAKE_COMMAND} .
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/language_tests/${language}
      RESULT_VARIABLE return_code
      OUTPUT_QUIET
      ERROR_QUIET
      )
    if(return_code EQUAL 0)
      set(${language_works} ON PARENT_SCOPE)
    else(return_code EQUAL 0)
      set(${language_works} OFF PARENT_SCOPE)
    endif(return_code EQUAL 0)
  else(return_code EQUAL 0)
    set(${language_works} OFF PARENT_SCOPE)
  endif(return_code EQUAL 0)
endfunction(workaround_9220)

# Temporary tests of the above function.
#workaround_9220(CXX CXX_language_works)
#message("CXX_language_works = ${CXX_language_works}")
#workaround_9220(CXXp CXXp_language_works)
#message("CXXp_language_works = ${CXXp_language_works}")
