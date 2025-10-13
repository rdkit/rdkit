add_executable(bench EXCLUDE_FROM_ALL smiles.cpp)
target_link_libraries(bench rdkitCatch GraphMol SmilesParse)

if(RDK_BUILD_CPP_TESTS)
  # add a fast version of the benchmarks to the default unit tests
  # to protect the benchmarks from bit-rot
  add_test(quick-bench bench --benchmark-samples 1 --benchmark-warmup-time 0)
endif(RDK_BUILD_CPP_TESTS)
