rdkit_library(RGroupDecomposition RGroupDecomp.cpp
  LINK_LIBRARIES
  FMCS ChemTransforms SubstructMatch SmilesParse 
  ${RDKit_THREAD_LIBS})

rdkit_headers(
   RGroupDecomp.h
   DEST GraphMol/RGroupDecomposition)


add_subdirectory(Wrap)


rdkit_test(testRGroupDecomp testRGroupDecomp.cpp 
           LINK_LIBRARIES RGroupDecomposition ChemTransforms FMCS FileParsers
           SmilesParse SubstructMatch GraphMol RDGeometryLib RDGeneral
           ${RDKit_THREAD_LIBS})
