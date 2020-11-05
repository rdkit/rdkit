
rdkit_library(RGroupDecomposition RGroupDecomp.cpp RGroupUtils.cpp
				  RGroupDecompParams.cpp RGroupScore.cpp RGroupFingerprintScore.cpp RGroupGa.cpp
  LINK_LIBRARIES
  FMCS ChemTransforms SubstructMatch SmilesParse Fingerprints
 GraphMol RDGeneral ga)
target_compile_definitions(RGroupDecomposition PRIVATE RDKIT_RGROUPDECOMPOSITION_BUILD)

rdkit_headers(
   RGroupDecomp.h
   DEST GraphMol/RGroupDecomposition)


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_test(testRGroupDecomp testRGroupDecomp.cpp
		LINK_LIBRARIES RGroupDecomposition )

add_executable(gaExample GaExample.cpp)
target_link_libraries(gaExample RGroupDecomposition)