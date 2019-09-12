rdkit_python_extension(rdFingerprintGenerator
                       FingerprintGeneratorWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Fingerprints 
                       DataStructs RDGeneral RDBoost)

rdkit_python_extension(rdMHFPFingerprint
                       MHFPWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Fingerprints 
                       DataStructs RDGeneral RDBoost)
                       
add_pytest(pyTestGenerator ${CMAKE_CURRENT_SOURCE_DIR}/testGenerators.py)

