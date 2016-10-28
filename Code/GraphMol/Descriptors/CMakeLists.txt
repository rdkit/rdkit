if(RDK_BUILD_DESCRIPTORS3D)
  set(DESC3D_HDRS MolDescriptors3D.h PBF.h PMI.h RDF.h MORSE.h WHIM.h GETAWAY.h) # auto3Dcorr.cpp)
  set(DESC3D_SOURCES PBF.cpp PMI.cpp RDF.cpp MORSE.cpp WHIM.cpp GETAWAY.cpp) # auto3Dcorr.h)
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
 rdkit_test(testRDF testRDF.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testMORSE testMORSE.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(test3D test3D.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testWHIM testWHIM.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testGETAWAY testGETAWAY.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
#  rdkit_test(test3Dauto test3Dauto.cpp
#LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )

endif(RDK_BUILD_DESCRIPTORS3D)



add_subdirectory(Wrap)
