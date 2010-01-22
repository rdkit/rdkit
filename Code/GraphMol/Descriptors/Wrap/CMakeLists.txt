rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES Descriptors FileParsers SmilesParse Subgraphs SubstructMatch GraphMol DataStructs RDGeneral RDGeometryLib RDBoost Fingerprints  )

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)





