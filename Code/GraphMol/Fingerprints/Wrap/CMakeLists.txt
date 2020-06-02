rdkit_python_extension(rdFingerprintGenerator
                       FingerprintGeneratorWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Fingerprints )

rdkit_python_extension(rdMHFPFingerprint
                       MHFPWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Fingerprints )

add_pytest(pyTestGenerator ${CMAKE_CURRENT_SOURCE_DIR}/testGenerators.py)
add_pytest(pyTestMHFP ${CMAKE_CURRENT_SOURCE_DIR}/testMHFP.py)
