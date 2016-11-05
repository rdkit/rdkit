function(build_minimalexample MINIMALEXAMPLE_EXECUTABLE LIBRARY)
  add_executable("${MINIMALEXAMPLE_EXECUTABLE}" MinimalExample.c)
  set (EXTRA_LIBS "${LIBRARY}")
  set_property(TARGET "${MINIMALEXAMPLE_EXECUTABLE}" PROPERTY C_STANDARD 90)
  target_link_libraries ("${MINIMALEXAMPLE_EXECUTABLE}" ${EXTRA_LIBS})
endfunction()

build_minimalexample("${SHARED_MINIMALEXAMPLE}" "${SHARED_LIBRARY}")
build_minimalexample("${STATIC_MINIMALEXAMPLE}" "${STATIC_LIBRARY}")
