find_package (Threads)
rdkit_library(Fingerprints
              Fingerprints.cpp MorganFingerprints.cpp AtomPairs.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse GraphMol ${CMAKE_THREAD_LIBS_INIT} )

rdkit_headers(AtomPairs.h
              Fingerprints.h
              MorganFingerprints.h DEST GraphMol/Fingerprints)

rdkit_test(testFingerprints test1.cpp LINK_LIBRARIES 
           Fingerprints_static FileParsers SubstructMatch SmilesParse
           Subgraphs_static GraphMol DataStructs RDGeometryLib
           RDGeneral  )


