
rdkit_library(GeneralizedSubstruct 
              XQMol.cpp
              LINK_LIBRARIES TautomerQuery SubstructMatch SmilesParse GraphMol) 
target_compile_definitions(GeneralizedSubstruct PRIVATE RDKIT_GENERALIZEDSUBSTRUCT_BUILD)

rdkit_headers(XQMol.h DEST GraphMol/GeneralizedSubstruct)

rdkit_catch_test(generalizedSubstructCatch catch_tests.cpp LINK_LIBRARIES FileParsers GeneralizedSubstruct)
