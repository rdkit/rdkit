if (PYTHON AND BUILD_PYTHON_WRAPPER)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/py_rdl/__init__.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
  set(PYTHON_INCLUDE_DIR "${LIB_SOURCE_DIR}")
  set(PYTHON_LIBRARY_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

  file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" NATIVE_CURRENT_SOURCE_DIR)

  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  add_custom_command(OUTPUT "${OUTPUT}"
                     COMMAND "${PYTHON}" "${SETUP_PY}" build "$<TARGET_FILE:RingDecomposerLibStatic>"
                     COMMAND ${CMAKE_COMMAND} -E touch "${OUTPUT}"
                     DEPENDS ${DEPS})

  add_custom_target(PYTHON_BUILD_TARGET ALL DEPENDS "${OUTPUT}")

  add_custom_command(TARGET PYTHON_BUILD_TARGET POST_BUILD
                     COMMAND "${PYTHON}" "${SETUP_PY}" install ${PYTHON_FLAGS} "$<TARGET_FILE:RingDecomposerLibStatic>")
endif()
