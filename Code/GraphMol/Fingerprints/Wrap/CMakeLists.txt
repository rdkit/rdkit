rdkit_python_extension(rdAtomPairGenerator
                       AtomPairGeneratorWrapper.cpp
                       DEST Chem
                       LINK_LIBRARIES Descriptors PartialCharges Fingerprints FileParsers SmilesParse Subgraphs
SubstructMatch MolTransforms GraphMol EigenSolvers DataStructs RDGeneral RDGeometryLib RDBoost)
                       
add_pytest(pyAtomPairGenerator ${CMAKE_CURRENT_SOURCE_DIR}/testAtomPair.py)

