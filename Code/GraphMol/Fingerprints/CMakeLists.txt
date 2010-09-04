rdkit_library(Fingerprints
              Fingerprints.cpp MorganFingerprints.cpp AtomPairs.cpp
              LINK_LIBRARIES Subgraphs)

rdkit_headers(AtomPairs.h
              Fingerprints.h
              MorganFingerprints.h DEST GraphMol/Fingerprints)

rdkit_test(testFingerprints test1.cpp LINK_LIBRARIES Fingerprints FileParsers SmilesParse Subgraphs GraphMol DataStructs RDGeometryLib RDGeneral  )

