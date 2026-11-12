
if(RDK_USE_HIGHWAY)
  set(BFGS_KERNELS BFGSKernelsHwy.cpp)
  # foreach_target.h re-includes the file, so it can't be merged with others
  set_source_files_properties(BFGSKernelsHwy.cpp PROPERTIES
                              SKIP_UNITY_BUILD_INCLUSION ON)
else()
  set(BFGS_KERNELS BFGSKernels.cpp)
endif()

rdkit_library(Optimizer
              BFGSOpt.cpp ${BFGS_KERNELS} LinearSearch.cpp
              LINK_LIBRARIES RDGeometryLib Trajectory RDGeneral)
target_compile_definitions(Optimizer PRIVATE RDKIT_OPTIMIZER_BUILD)
if(RDK_USE_HIGHWAY)
  target_link_libraries(Optimizer PRIVATE hwy::hwy)
  if(TARGET Optimizer_static)
    target_link_libraries(Optimizer_static PRIVATE hwy::hwy)
  endif()
endif()

rdkit_headers(BFGSOpt.h BFGSKernels.h DEST Numerics/Optimizer)

rdkit_catch_test(testOptimizer testOptimizer.cpp LINK_LIBRARIES Optimizer )
if(RDK_USE_HIGHWAY AND TARGET testOptimizer)
  # lets the test run the kernels on every SIMD target the CPU supports
  target_link_libraries(testOptimizer PRIVATE hwy::hwy)
  target_compile_definitions(testOptimizer PRIVATE RDK_USE_HIGHWAY)
endif()
