
rdkit_library(GenericGroups 
              GenericGroups.cpp 
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(GenericGroups PRIVATE RDKIT_GENERICGROUPS_BUILD)

rdkit_headers(GenericGroups.h DEST GraphMol/GenericGroups)

rdkit_catch_test(testGenericGroups generic_tests.cpp 
    LINK_LIBRARIES GenericGroups FileParsers SmilesParse SubstructMatch)
