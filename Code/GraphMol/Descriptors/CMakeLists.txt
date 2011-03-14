rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp
              LINK_LIBRARIES SmilesParse FileParsers Subgraphs SubstructMatch)

rdkit_headers(Crippen.h Lipinski.h
              MolDescriptors.h
              MolSurf.h DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp 
LINK_LIBRARIES Descriptors FileParsers SmilesParse Subgraphs SubstructMatch GraphMol DataStructs RDGeneral RDGeometryLib  )

add_subdirectory(Wrap)



