rdkit_python_extension(rdAtomPairGenerator
                       FingerprintGeneratorWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Fingerprints 
                       DataStructs RDGeneral RDBoost)
                       
add_pytest(pyAtomPairGenerator ${CMAKE_CURRENT_SOURCE_DIR}/testAtomPair.py)

