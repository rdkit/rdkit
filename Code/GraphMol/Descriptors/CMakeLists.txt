rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp AtomPairs.cpp MolSurf.cpp
              LINK_LIBRARIES FileParsers Subgraphs SubstructMatch)

rdkit_test(testDescriptors test.cpp LINK_LIBRARIES Descriptors)

add_subdirectory(Wrap)



