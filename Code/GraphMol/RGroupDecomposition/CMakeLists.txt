
rdkit_library(RGroupDecomposition RGroupDecomp.cpp RGroupDecompData.cpp
  RGroupData.cpp RGroupUtils.cpp RGroupCore.cpp
  RGroupDecompParams.cpp RGroupScore.cpp RGroupFingerprintScore.cpp
  RGroupGa.cpp RGroupDecompJSONParsers.cpp
  LINK_LIBRARIES
  FMCS ChemTransforms SubstructMatch SmilesParse Fingerprints
    GraphMol RDGeneral ga MolEnumerator TautomerQuery)
target_compile_definitions(RGroupDecomposition PRIVATE RDKIT_RGROUPDECOMPOSITION_BUILD)

rdkit_headers(
  RGroupDecomp.h RGroupDecompParams.h RGroupDecompJSONParsers.h
  DEST GraphMol/RGroupDecomposition)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

rdkit_catch_test(testRGroupDecomp testRGroupDecomp.cpp
  LINK_LIBRARIES RGroupDecomposition)

rdkit_catch_test(testRGroupDecompInternals testRGroupInternals.cpp
  LINK_LIBRARIES RGroupDecomposition)

rdkit_catch_test(rgroupCatchTests catch_rgd.cpp 
  LINK_LIBRARIES RGroupDecomposition)

if (RDK_BUILD_CPP_TESTS AND NOT RDK_BUILD_CFFI_LIB AND NOT RDK_BUILD_MINIMAL_LIB)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS program_options CONFIG)
  add_executable(gaExample GaExample.cpp)
  if(NOT Boost_USE_STATIC_LIBS)
    target_compile_definitions(gaExample PUBLIC -DBOOST_PROGRAM_OPTIONS_DYN_LINK)
  endif()
  target_link_libraries(gaExample RGroupDecomposition ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()
