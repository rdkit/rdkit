
rdkit_library(ChemTransforms ChemTransforms.cpp
  MolFragmenter.cpp MolFragmenterJSONParser.cpp LINK_LIBRARIES
  GraphMol SubstructMatch SmilesParse Depictor)
target_compile_definitions(ChemTransforms PRIVATE -DRDKIT_CHEMTRANSFORMS_BUILD)

rdkit_headers(ChemTransforms.h
   MolFragmenter.h MolFragmenterJSONParser.h
   DEST GraphMol/ChemTransforms)

# there's no Wrap subdirectory on the main trunk (but in "minimal" there is)..
#add_subdirectory(Wrap)

rdkit_test(testChemTransforms testChemTransforms.cpp 
           LINK_LIBRARIES ChemTransforms FileParsers)

rdkit_catch_test(chemTransformsTestCatch catch_tests.cpp LINK_LIBRARIES ChemTransforms FileParsers)
