if(RDK_BUILD_DESCRIPTORS3D)

  set(DESC3D_HDRS MolDescriptors3D.h EEM.h PBF.h PMI.h AUTOCORR3D.h RDF.h MORSE.h GETAWAY.h WHIM.h CoulombMat.h)
  set(DESC3D_SOURCES EEM.cpp PBF.cpp PMI.cpp AUTOCORR3D.cpp RDF.cpp MORSE.cpp GETAWAY.cpp WHIM.cpp CoulombMat.cpp)

endif(RDK_BUILD_DESCRIPTORS3D)

if(RDK_BUILD_OSMORDRED_SUPPORT)
 find_package(LAPACK REQUIRED)
 find_package(LAPACKE REQUIRED)

 if(NOT LAPACKE_FOUND)
   message(fatal_error "-- LAPACKE not found")
 endif()

 include_directories(Descriptors LAPACK::LAPACK ${LAPACK_INCLUDE_DIRS})

endif(RDK_BUILD_OSMORDRED_SUPPORT)

rdkit_library(Descriptors
              Crippen.cpp BCUT.cpp MolDescriptors.cpp MolSurf.cpp Lipinski.cpp ConnectivityDescriptors.cpp
              MQN.cpp
              Property.cpp
              AUTOCORR2D.cpp Data3Ddescriptors.cpp MolData3Ddescriptors.cpp
              USRDescriptor.cpp AtomFeat.cpp
              OxidationNumbers.cpp
              DCLV.cpp
	      Osmordred.cpp
              ${DESC3D_SOURCES}
              LINK_LIBRARIES DataStructs Fingerprints PartialCharges SmilesParse FileParsers Subgraphs SubstructMatch MolTransforms GraphMol
                 EigenSolvers RDGeneral)
target_compile_definitions(Descriptors PRIVATE RDKIT_DESCRIPTORS_BUILD)
target_link_libraries(Descriptors PUBLIC ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES} ${BLAS_LIBRARIES})
rdkit_headers(Crippen.h BCUT.h Lipinski.h
              MolDescriptors.h
              MolSurf.h
              ConnectivityDescriptors.h MQN.h
              RegisterDescriptor.h
              Property.h
              AUTOCORR2D.h Data3Ddescriptors.h MolData3Ddescriptors.h
              USRDescriptor.h AtomFeat.h
              OxidationNumbers.h
              DCLV.h
	      Osmordred.h
              ${DESC3D_HDRS}
              DEST GraphMol/Descriptors)

rdkit_test(testDescriptors test.cpp
LINK_LIBRARIES PartialCharges Descriptors )

rdkit_test(testAUTOCORR2D testAUTOCORR2D.cpp
LINK_LIBRARIES  Descriptors )

rdkit_test(testATF testATF.cpp
LINK_LIBRARIES  Descriptors )

if(RDK_HAS_EIGEN)
  rdkit_test(testBCUT testBCUT.cpp
  LINK_LIBRARIES  Descriptors )
endif(RDK_HAS_EIGEN)

if(RDK_BUILD_DESCRIPTORS3D)
  rdkit_test(testEEM testEEM.cpp
LINK_LIBRARIES  Descriptors )
  rdkit_test(testPBF testPBF.cpp
LINK_LIBRARIES  Descriptors )
 rdkit_test(testRDFcustom testRDFcustom.cpp
LINK_LIBRARIES  Descriptors )
 rdkit_test(testRDF testRDF.cpp
LINK_LIBRARIES  Descriptors )
 rdkit_test(testMORSE testMORSE.cpp
LINK_LIBRARIES  Descriptors )
rdkit_test(test3D test3D.cpp
LINK_LIBRARIES  Descriptors )
  rdkit_test(testWHIM testWHIM.cpp
LINK_LIBRARIES  Descriptors )
  rdkit_test(testGETAWAY testGETAWAY.cpp
LINK_LIBRARIES  Descriptors )
  rdkit_test(testAUTOCORR3D testAUTOCORR3D.cpp
LINK_LIBRARIES  Descriptors )
 rdkit_test(testCoulombMat testCoulombMat.cpp
LINK_LIBRARIES  Descriptors ForceFieldHelpers DistGeomHelpers  )
endif(RDK_BUILD_DESCRIPTORS3D)

rdkit_catch_test(descriptorsTestCatch catch_tests.cpp LINK_LIBRARIES Descriptors SmilesParse FileParsers)

rdkit_catch_test(testOsmordred test_osmordred.cpp LINK_LIBRARIES Descriptors SmilesParse FileParsers)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
