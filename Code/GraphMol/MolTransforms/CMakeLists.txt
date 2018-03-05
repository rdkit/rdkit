rdkit_library(MolTransforms MolTransforms.cpp
              LINK_LIBRARIES GraphMol EigenSolvers)

rdkit_headers(MolTransforms.h DEST GraphMol/MolTransforms)

rdkit_test(testMolTransforms test1.cpp 
           LINK_LIBRARIES MolTransforms 
FileParsers SmilesParse 
EigenSolvers 
GraphMol RDGeometryLib RDGeneral 
)

add_subdirectory(Wrap)
