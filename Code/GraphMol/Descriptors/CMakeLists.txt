rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp
              LINK_LIBRARIES FileParsers Subgraphs SubstructMatch)

rdkit_headers(Crippen.h
              MolDescriptors.h
              MolSurf.h DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp 
LINK_LIBRARIES Descriptors FileParsers SmilesParse Subgraphs SubstructMatch GraphMol DataStructs RDGeneral RDGeometryLib  )

add_subdirectory(Wrap)



