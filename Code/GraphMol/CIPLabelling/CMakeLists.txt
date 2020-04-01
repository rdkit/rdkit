
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_CIPLABELLING_BUILD)
rdkit_library(CIPLabelling
              RDKitCIPMol.cpp RDKitCIPLabeller.cpp
              LINK_LIBRARIES GraphMol)

rdkit_headers(CIPLabelling.h
              DEST GraphMol/CIPLabelling)

rdkit_catch_test(testCIPLabelling catch_tests.cpp RDKitCIPMol.cpp
                 LINK_LIBRARIES CIPLabelling SmilesParse GraphMol)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()









