include(BoostUtils.cmake)
macro(rdkit_python_extension)
  PARSE_ARGUMENTS(RDKPY
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN}
    )
  CAR(RDKPY_NAME ${RDKPY_DEFAULT_ARGS})
  CDR(RDKPY_SOURCES ${RDKPY_DEFAULT_ARGS})
  PYTHON_ADD_MODULE(${RDKPY_NAME} ${RDKPY_SOURCES}  )
  set_target_properties(${RDKPY_NAME} PROPERTIES PREFIX "")
if(MSVC)
  set_target_properties(${RDKPY_NAME} PROPERTIES SUFFIX ".pyd")
endif(MSVC)  
  target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

  INSTALL(TARGETS ${RDKPY_NAME} 
      LIBRARY DESTINATION ${RDKit_PythonDir}/${RDKPY_DEST} )
endmacro(rdkit_python_extension)

macro(rdkit_test)
  PARSE_ARGUMENTS(RDKTEST
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN}
    )
  CAR(RDKTEST_NAME ${RDKTEST_DEFAULT_ARGS})
  CDR(RDKTEST_SOURCES ${RDKTEST_DEFAULT_ARGS})
  add_executable(${RDKTEST_NAME} ${RDKTEST_SOURCES})
  target_link_libraries(${RDKTEST_NAME} ${RDKTEST_LINK_LIBRARIES} )
  add_test(${RDKTEST_NAME} ${EXECUTABLE_OUTPUT_PATH}/${RDKTEST_NAME} )
endmacro(rdkit_test)

