rdkit_library(ChemTransforms ChemTransforms.cpp MolFragmenter.cpp LINK_LIBRARIES
SubstructMatch SmilesParse )

rdkit_headers(ChemTransforms.h DEST GraphMol/ChemTransforms)

# there's no Wrap subdirectory on the main trunk (but in "minimal" there is)..
#add_subdirectory(Wrap)

rdkit_test(testChemTransforms testChemTransforms.cpp 
           LINK_LIBRARIES ChemTransforms FileParsers SmilesParse SubstructMatch GraphMol RDGeometryLib RDGeneral )



