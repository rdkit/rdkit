include(BoostUtils)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac OS X specific code
  set(RDKit_VERSION "${RDKit_Year}.${RDKit_Month}")
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(RDKit_VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(RDKit_RELEASENAME "${RDKit_Year}.${RDKit_Month}")
if (RDKit_Revision)
  set(RDKit_RELEASENAME "${RDKit_RELEASENAME}.${RDKit_Revision}")
  set(RDKit_VERSION "${RDKit_VERSION}.${RDKit_Revision}")
else(RDKit_Revision)
  set(RDKit_VERSION "${RDKit_VERSION}.0")
endif(RDKit_Revision)

set(RDKit_EXPORTED_TARGETS rdkit-targets)

macro(rdkit_library)
  PARSE_ARGUMENTS(RDKLIB
    "LINK_LIBRARIES;DEST"
    "SHARED"
    ${ARGN})
  CAR(RDKLIB_NAME ${RDKLIB_DEFAULT_ARGS})
  CDR(RDKLIB_SOURCES ${RDKLIB_DEFAULT_ARGS})
  if(MSVC)
    add_library(${RDKLIB_NAME} ${RDKLIB_SOURCES})
    target_link_libraries(${RDKLIB_NAME} ${Boost_SYSTEM_LIBRARY} )
    INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
            DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
            COMPONENT dev )
  else(MSVC)
    # we're going to always build in shared mode since we
    # need exceptions to be (correctly) catchable across
    # boundaries. As of now (June 2010), this doesn't work
    # with g++ unless libraries are shared.
      add_library(${RDKLIB_NAME} SHARED ${RDKLIB_SOURCES})
      INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
              DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
              COMPONENT runtime )
      if(RDK_INSTALL_STATIC_LIBS)
        add_library(${RDKLIB_NAME}_static ${RDKLIB_SOURCES})
        INSTALL(TARGETS ${RDKLIB_NAME}_static EXPORT ${RDKit_EXPORTED_TARGETS}
                DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
                COMPONENT dev )
      endif(RDK_INSTALL_STATIC_LIBS)
    IF(RDKLIB_LINK_LIBRARIES)
      target_link_libraries(${RDKLIB_NAME} ${RDKLIB_LINK_LIBRARIES})
    ENDIF(RDKLIB_LINK_LIBRARIES)
  endif(MSVC)
  if(WIN32)
    set_target_properties(${RDKLIB_NAME} PROPERTIES 
                          OUTPUT_NAME "${RDKLIB_NAME}" 
                          VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}")
  else(WIN32)
    set_target_properties(${RDKLIB_NAME} PROPERTIES 
                          OUTPUT_NAME ${RDKLIB_NAME} 
                          VERSION ${RDKit_VERSION} 
                          SOVERSION ${RDKit_ABI} )
  endif(WIN32)			  
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
    install(FILES ${RDKHDR_DEFAULT_ARGS}
            DESTINATION ${RDKit_HdrDir}/${RDKHDR_DEST}
            COMPONENT dev )
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
                          ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} )

    INSTALL(TARGETS ${RDKPY_NAME} 
            LIBRARY DESTINATION ${RDKit_PythonDir}/${RDKPY_DEST})
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
  endif(RDK_BUILD_PYTHON_WRAPPERS)
endmacro(add_pytest)
