if(RDK_BUILD_DESCRIPTORS3D)
  set(DESC3D_HDRS MolDescriptors3D.h PBF.h PMI.h)
  set(DESC3D_SOURCES PBF.cpp PMI.cpp)
endif(RDK_BUILD_DESCRIPTORS3D)

rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp ConnectivityDescriptors.cpp
              MQN.cpp
              Property.cpp
              ${DESC3D_SOURCES}
              LINK_LIBRARIES PartialCharges SmilesParse FileParsers Subgraphs SubstructMatch MolTransforms GraphMol
                 EigenSolvers ${RDKit_THREAD_LIBS})



rdkit_headers(Crippen.h Lipinski.h
              MolDescriptors.h
              MolSurf.h
              ConnectivityDescriptors.h MQN.h
              RegisterDescriptor.h
              Property.h
              ${DESC3D_HDRS}
              DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp
LINK_LIBRARIES PartialCharges Descriptors FileParsers SmilesParse Subgraphs SubstructMatch MolTransforms GraphMol EigenSolvers DataStructs RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )

if(RDK_BUILD_DESCRIPTORS3D)
  rdkit_test(testPBF testPBF.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
rdkit_test(test3D test3D.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
endif(RDK_BUILD_DESCRIPTORS3D)



add_subdirectory(Wrap)
