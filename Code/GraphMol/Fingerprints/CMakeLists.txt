rdkit_library(Fingerprints
              Fingerprints.cpp PatternFingerprints.cpp MorganFingerprints.cpp AtomPairs.cpp MACCS.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse GraphMol
                ${RDKit_THREAD_LIBS} )

rdkit_headers(AtomPairs.h
              Fingerprints.h
              MorganFingerprints.h
              MACCS.h
              DEST GraphMol/Fingerprints)

rdkit_test(testFingerprints test1.cpp LINK_LIBRARIES 
           Fingerprints FileParsers SubstructMatch SmilesParse
           Subgraphs GraphMol DataStructs RDGeometryLib
           RDGeneral  ${RDKit_THREAD_LIBS}  )


