
rdkit_library(ChemTransforms ChemTransforms.cpp
  MolFragmenter.cpp MolFragmenterJSONParser.cpp MolZip.cpp LINK_LIBRARIES
  GraphMol SubstructMatch SmilesParse Depictor)
target_compile_definitions(ChemTransforms PRIVATE -DRDKIT_CHEMTRANSFORMS_BUILD)

# MolFragmenterJSONParser.cpp #defines USE_BETTER_ENUMS and so is incompatible with
# the plain-enum precompiled header; compile it without the PCH (see BetterEnums.h).
if(RDK_USE_PRECOMPILED_HEADERS)
  set_source_files_properties(MolFragmenterJSONParser.cpp PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
endif()

rdkit_headers(ChemTransforms.h
   MolFragmenter.h MolFragmenterJSONParser.h
   DEST GraphMol/ChemTransforms)

# there's no Wrap subdirectory on the main trunk (but in "minimal" there is)..
#add_subdirectory(Wrap)

rdkit_test(testChemTransforms testChemTransforms.cpp 
           LINK_LIBRARIES ChemTransforms FileParsers)

rdkit_catch_test(chemTransformsTestCatch catch_tests.cpp LINK_LIBRARIES ChemTransforms FileParsers)
