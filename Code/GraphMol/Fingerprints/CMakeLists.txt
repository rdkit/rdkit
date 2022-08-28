
rdkit_library(Fingerprints
              Fingerprints.cpp PatternFingerprints.cpp MorganFingerprints.cpp
              AtomPairs.cpp MACCS.cpp MHFP.cpp FingerprintGenerator.cpp 
              AtomPairGenerator.cpp MorganGenerator.cpp RDKitFPGenerator.cpp 
              FingerprintUtil.cpp TopologicalTorsionGenerator.cpp
              LINK_LIBRARIES DataStructs Subgraphs SubstructMatch SmilesParse GraphMol RDGeneral
              )
target_compile_definitions(Fingerprints PRIVATE RDKIT_FINGERPRINTS_BUILD)

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
           Fingerprints FileParsers )

rdkit_test(testMHFPFingerprint testMHFPFingerprint.cpp LINK_LIBRARIES 
           Fingerprints FileParsers )

rdkit_test(testFingerprintGenerators testFingerprintGenerators.cpp LINK_LIBRARIES 
           Fingerprints FileParsers )

rdkit_catch_test(fpTestCatch catch_tests.cpp LINK_LIBRARIES Fingerprints FileParsers)


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()