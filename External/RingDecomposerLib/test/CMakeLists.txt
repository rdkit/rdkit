file(GLOB files "${TEST_DIR}/*.dimacs")
foreach(filename ${files})
    get_filename_component(tname "${filename}" NAME_WE)
    add_test(NAME "test_${tname}" WORKING_DIRECTORY "${TEST_DIR}"
             COMMAND "${STATIC_TEST}" "validate" "${filename}" "200")
endforeach()

add_test(NAME "test_minimal_static" WORKING_DIRECTORY "${TEST_DIR}"
         COMMAND "${STATIC_MINIMALEXAMPLE}")

if (PYTHON AND BUILD_PYTHON_WRAPPER)
    foreach(filename ${files})
        get_filename_component(tname "${filename}" NAME_WE)
        add_test(NAME "python_test_${tname}" WORKING_DIRECTORY "${TEST_DIR}"
	             COMMAND "${PYTHON}" "${PYTHON_TEST}" "--input" "${filename}")
    endforeach()
endif()
