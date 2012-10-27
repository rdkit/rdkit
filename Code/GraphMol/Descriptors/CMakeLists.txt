find_package (Threads)

rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp ConnectivityDescriptors.cpp
              LINK_LIBRARIES PartialCharges SmilesParse FileParsers Subgraphs SubstructMatch 
                ${Boost_THREAD_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

rdkit_headers(Crippen.h Lipinski.h
              MolDescriptors.h
              MolSurf.h DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp 
LINK_LIBRARIES PartialCharges Descriptors FileParsers SmilesParse Subgraphs SubstructMatch GraphMol DataStructs RDGeneral RDGeometryLib )


add_subdirectory(Wrap)



