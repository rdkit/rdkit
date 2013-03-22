rdkit_library(ChemTransforms ChemTransforms.cpp MolFragmenter.cpp LINK_LIBRARIES
  SubstructMatch SmilesParse 
  ${RDKit_THREAD_LIBS})

rdkit_headers(ChemTransforms.h
   MolFragmenter.h
   DEST GraphMol/ChemTransforms)

# there's no Wrap subdirectory on the main trunk (but in "minimal" there is)..
#add_subdirectory(Wrap)

#rdkit_test(testChemTransforms testChemTransforms.cpp 
#           LINK_LIBRARIES ChemTransforms FileParsers SmilesParse SubstructMatch GraphMol RDGeometryLib RDGeneral )
rdkit_test(testChemTransforms testChemTransforms.cpp 
           LINK_LIBRARIES ChemTransforms FileParsers SmilesParse SubstructMatch
           GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS})



