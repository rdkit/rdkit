remove_definitions(-DRDKIT_DESCRIPTORS_BUILD)
rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES Descriptors PartialCharges Fingerprints FileParsers SmilesParse Subgraphs
SubstructMatch MolTransforms GraphMol EigenSolvers DataStructs RDGeneral RDGeometryLib RDBoost )

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)
add_pytest(pyMolDescriptors3D ${CMAKE_CURRENT_SOURCE_DIR}/test3D.py)
