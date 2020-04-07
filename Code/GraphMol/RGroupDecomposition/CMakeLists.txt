
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_RGROUPDECOMPOSITION_BUILD)
rdkit_library(RGroupDecomposition RGroupDecomp.cpp
  LINK_LIBRARIES
  FMCS ChemTransforms SubstructMatch SmilesParse 
 GraphMol RDGeneral)

rdkit_headers(
   RGroupDecomp.h
   DEST GraphMol/RGroupDecomposition)


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()


rdkit_test(testRGroupDecomp testRGroupDecomp.cpp 
           LINK_LIBRARIES RGroupDecomposition )
