add_executable(bench EXCLUDE_FROM_ALL bench_common.cpp
  descriptors.cpp
  fingerprint.cpp
  meta.cpp
  mol.cpp
  molops.cpp
  pickle.cpp
  smiles.cpp
  stereo.cpp
  substruct_match.cpp
)
target_link_libraries(bench rdkitCatch
  CIPLabeler
  Descriptors
  Fingerprints
  SmilesParse
)

if(RDK_BUILD_INCHI_SUPPORT)
  target_link_libraries(bench RDInchiLib)
  target_sources(bench PRIVATE inchi.cpp)
endif(RDK_BUILD_INCHI_SUPPORT)

if(RDK_BUILD_CPP_TESTS)
  # add a fast version of the benchmarks to the default unit tests
  # to protect the benchmarks from bit-rot
  add_test(
    NAME quickbench
    COMMAND bench --benchmark-samples 1 --benchmark-warmup-time 0
  )
  
  # ctest does not automatically build test executables
  # so, like with other unit tests in this repo, we need to manually
  # build the target before running ctest. This usually gets done with
  # an `all` or an `install` build. We re-add bench to he `all` group
  # to keep it similar to the other unit tests
  set_target_properties(bench PROPERTIES EXCLUDE_FROM_ALL FALSE)
endif(RDK_BUILD_CPP_TESTS)
