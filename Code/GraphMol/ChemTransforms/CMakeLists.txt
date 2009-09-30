add_library(ChemTransforms ChemTransforms.cpp)

add_subdirectory(Wrap)

rdkit_test(testChemTransforms testChemTransforms.cpp LINK_LIBRARIES
  ChemTransforms SubstructMatch SmilesParse FileParsers GraphMol RDGeneral RDGeometryLib )



