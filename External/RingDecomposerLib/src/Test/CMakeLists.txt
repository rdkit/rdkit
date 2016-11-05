function(build_test TEST_EXECUTABLE LIBRARY)
  add_executable("${TEST_EXECUTABLE}" Test.c
                 TestDemo.c TestValidate.c)
  set (EXTRA_LIBS "${LIBRARY}")
  set_property(TARGET "${TEST_EXECUTABLE}" PROPERTY C_STANDARD 90)
  target_link_libraries ("${TEST_EXECUTABLE}" ${EXTRA_LIBS})
endfunction()

build_test("${SHARED_TEST}" "${SHARED_LIBRARY}")
build_test("${STATIC_TEST}" "${STATIC_LIBRARY}")
