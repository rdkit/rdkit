rdkit_library(SubstructMatch 
              SubstructMatch.cpp SubstructUtils.cpp
              LINK_LIBRARIES GraphMol)

rdkit_test(testSubstructMatch test1.cpp LINK_LIBRARIES SubstructMatch
GraphMol RDGeometryLib RDGeneral )

