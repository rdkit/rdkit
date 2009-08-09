add_library(SubstructMatch SubstructMatch.cpp SubstructUtils.cpp)

rdkit_test(substructMatchTest1 test1.cpp LINK_LIBRARIES SubstructMatch GraphMol RDGeneral RDGeometry)


