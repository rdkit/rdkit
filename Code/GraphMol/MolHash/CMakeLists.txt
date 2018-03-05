rdkit_library(MolHash
              MolHash.cpp HashToString.cpp
              LINK_LIBRARIES Depictor Descriptors GraphMol RDGeneral)

rdkit_headers(MolHash.h             
              DEST GraphMol/MolHash)

rdkit_test(testMolHash testMolHash.cpp 
LINK_LIBRARIES MolHash GraphMol ChemTransforms Depictor Descriptors FileParsers SmilesParse RDGeneral RDGeometryLib)

add_subdirectory(Wrap)

