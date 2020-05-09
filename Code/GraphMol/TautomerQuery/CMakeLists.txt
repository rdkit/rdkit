
# remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
# add_definitions(-DRDKIT_SUBSTRUCTMATCH_BUILD)
rdkit_library(tautomerQuery 
              TautomerQuery.cpp
              LINK_LIBRARIES GraphMol RDGeneral MolStandardize)

rdkit_headers(TautomerQuery.h
              DEST GraphMol/TautomerQuery)

rdkit_catch_test(tautomerQueryTestCatch catch_tests.cpp LINK_LIBRARIES tautomerQuery)
