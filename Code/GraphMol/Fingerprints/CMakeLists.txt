rdkit_library(Fingerprints
              Fingerprints.cpp MorganFingerprints.cpp
              LINK_LIBRARIES Subgraphs)

rdkit_test(testFingerprints test1.cpp LINK_LIBRARIES SmilesParse Subgraphs Fingerprints GraphMol DataStructs RDGeometryLib RDGeneral ${LAPACK_LIBRARIES} )

