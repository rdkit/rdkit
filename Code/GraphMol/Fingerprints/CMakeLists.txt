
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_FINGERPRINTS_BUILD)
rdkit_library(Fingerprints
              Fingerprints.cpp PatternFingerprints.cpp MorganFingerprints.cpp
              AtomPairs.cpp MACCS.cpp MHFP.cpp FingerprintGenerator.cpp 
              AtomPairGenerator.cpp MorganGenerator.cpp RDKitFPGenerator.cpp 
              FingerprintUtil.cpp TopologicalTorsionGenerator.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse GraphMol
              ${RDKit_THREAD_LIBS} )

rdkit_headers(AtomPairs.h
              Fingerprints.h
              MorganFingerprints.h
              MACCS.h
              MHFP.h
              FingerprintGenerator.h
              AtomPairGenerator.h
              MorganGenerator.h
              RDKitFPGenerator.h
              TopologicalTorsionGenerator.h
              FingerprintUtil.h
              DEST GraphMol/Fingerprints)

rdkit_test(testFingerprints test1.cpp LINK_LIBRARIES 
           Fingerprints FileParsers SubstructMatch SmilesParse
           Subgraphs GraphMol DataStructs RDGeometryLib
           RDGeneral  ${RDKit_THREAD_LIBS}  )

rdkit_test(testMHFPFingerprint testMHFPFingerprint.cpp LINK_LIBRARIES 
           Fingerprints SmilesParse GraphMol DataStructs FileParsers
           RDGeneral )

rdkit_test(testFingerprintGenerators testFingerprintGenerators.cpp LINK_LIBRARIES 
           Fingerprints SmilesParse GraphMol DataStructs FileParsers
           RDGeneral )
rdkit_catch_test(fpTestCatch catch_tests.cpp LINK_LIBRARIES Fingerprints SmilesParse GraphMol DataStructs FileParsers
           RDGeneral)


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()