include(CheckCXXSourceCompiles)

# Probe whether the host compiler accepts the AVX-512 VPOPCNTDQ flag set we
# need for the runtime-dispatched fast path. When it does, we build an extra
# translation unit that holds the AVX-512 kernels and a CPUID check at
# startup selects between scalar and AVX-512 at first use.
set(BULKSIM_SOURCES BulkSimilarityCpu.cpp)

if(MSVC)
  set(_bulksim_avx512_flag "/arch:AVX512")
else()
  set(_bulksim_avx512_flag "-mavx512f;-mavx512vpopcntdq")
endif()

set(_bulksim_save_flags ${CMAKE_REQUIRED_FLAGS})
if(MSVC)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} /arch:AVX512")
else()
  set(CMAKE_REQUIRED_FLAGS
      "${CMAKE_REQUIRED_FLAGS} -mavx512f -mavx512vpopcntdq")
endif()
check_cxx_source_compiles("
#include <immintrin.h>
int main() {
  __m512i v = _mm512_setzero_si512();
  __m512i p = _mm512_popcnt_epi64(v);
  return static_cast<int>(_mm512_reduce_add_epi64(p));
}" RDK_BULKSIM_COMPILER_HAS_AVX512_VPOPCNTQ)
set(CMAKE_REQUIRED_FLAGS ${_bulksim_save_flags})

if(RDK_BULKSIM_COMPILER_HAS_AVX512_VPOPCNTQ)
  list(APPEND BULKSIM_SOURCES BulkSimilarityCpuAvx512.cpp)
endif()

rdkit_library(BulkSimilarity
              ${BULKSIM_SOURCES}
              LINK_LIBRARIES DataStructs RDGeneral)
target_compile_definitions(BulkSimilarity PRIVATE RDKIT_BULKSIMILARITY_BUILD)

if(RDK_BULKSIM_COMPILER_HAS_AVX512_VPOPCNTQ)
  target_compile_definitions(BulkSimilarity
                             PRIVATE RDK_BULKSIM_HAS_AVX512_VPOPCNTQ)
  set_source_files_properties(BulkSimilarityCpuAvx512.cpp
                              PROPERTIES COMPILE_OPTIONS
                              "${_bulksim_avx512_flag}")
endif()

rdkit_headers(BulkSimilarity.h DEST DataStructs/BulkSimilarity)

rdkit_catch_test(testBulkSimilarity catch_tests.cpp
                 LINK_LIBRARIES BulkSimilarity DataStructs)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()
