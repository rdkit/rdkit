rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp ConnectivityDescriptors.cpp
              MQN.cpp 
              LINK_LIBRARIES PartialCharges SmilesParse FileParsers Subgraphs SubstructMatch 
                ${RDKit_THREAD_LIBS})

rdkit_headers(Crippen.h Lipinski.h
              MolDescriptors.h
              MolSurf.h
              ConnectivityDescriptors.h MQN.h
              DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp 
LINK_LIBRARIES PartialCharges Descriptors FileParsers SmilesParse Subgraphs SubstructMatch GraphMol DataStructs RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )


add_subdirectory(Wrap)



