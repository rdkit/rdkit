rdkit_library(DetermineBonds DetermineBonds.cpp
              LINK_LIBRARIES GenericGroups GraphMol RDGeneral FileParsers EHTLib)

target_compile_definitions(DetermineBonds PRIVATE RDKIT_DETERMINEBONDS_BUILD)

rdkit_headers(DetermineBonds.h DEST GraphMol/DetermineBonds)

rdkit_catch_test(determineBondsCatchTest catch_tests.cpp 
           LINK_LIBRARIES DetermineBonds FileParsers SmilesParse DistGeomHelpers ForceFieldHelpers)