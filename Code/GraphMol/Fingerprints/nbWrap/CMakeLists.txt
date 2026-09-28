remove_definitions(-DRDKIT_FINGERPRINTS_BUILD)
rdkit_nanobind_extension(rdFingerprintGenerator
                         FingerprintGeneratorWrapper.cpp
                         DEST Chem
                         LINK_LIBRARIES Fingerprints)

rdkit_nanobind_extension(rdMHFPFingerprint
                         MHFPWrapper.cpp
                         DEST Chem
                         LINK_LIBRARIES Fingerprints)

add_pytest(pyTestGenerator ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testGenerators.py)
add_pytest(pyTestMHFP ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMHFP.py)
