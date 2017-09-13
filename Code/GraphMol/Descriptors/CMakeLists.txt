if(RDK_BUILD_DESCRIPTORS3D)

  set(DESC3D_HDRS MolDescriptors3D.h PBF.h PMI.h Data3Ddescriptors.h MolData3Ddescriptors.h AUTOCORR3D.h AUTOCORR2D.h RDF.h MORSE.h GETAWAY.h WHIM.h)
  set(DESC3D_SOURCES PBF.cpp PMI.cpp Data3Ddescriptors.cpp MolData3Ddescriptors.cpp AUTOCORR3D.cpp AUTOCORR2D.cpp RDF.cpp MORSE.cpp GETAWAY.cpp WHIM.cpp)

endif(RDK_BUILD_DESCRIPTORS3D)

rdkit_library(Descriptors
              Crippen.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp ConnectivityDescriptors.cpp
              MQN.cpp
              Property.cpp
              USRDescriptor.cpp
              ${DESC3D_SOURCES}
              LINK_LIBRARIES PartialCharges SmilesParse FileParsers Subgraphs SubstructMatch MolTransforms GraphMol
                 EigenSolvers ${RDKit_THREAD_LIBS})

rdkit_headers(Crippen.h Lipinski.h
              MolDescriptors.h
              MolSurf.h
              ConnectivityDescriptors.h MQN.h
              RegisterDescriptor.h
              Property.h
              USRDescriptor.h
              ${DESC3D_HDRS}
              DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp
LINK_LIBRARIES PartialCharges Descriptors FileParsers SmilesParse Subgraphs SubstructMatch MolTransforms GraphMol EigenSolvers DataStructs RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )

if(RDK_BUILD_DESCRIPTORS3D)
  rdkit_test(testPBF testPBF.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
 rdkit_test(testRDF testRDF.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
 rdkit_test(testMORSE testMORSE.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
rdkit_test(test3D test3D.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testWHIM testWHIM.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testGETAWAY testGETAWAY.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testAUTOCORR3D testAUTOCORR3D.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )
  rdkit_test(testAUTOCORR2D testAUTOCORR2D.cpp
LINK_LIBRARIES  Descriptors FileParsers SmilesParse MolTransforms PartialCharges GraphMol DataStructs EigenSolvers RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS} )

endif(RDK_BUILD_DESCRIPTORS3D)


add_subdirectory(Wrap)
