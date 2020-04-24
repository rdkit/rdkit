
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_CIPLABELER_BUILD)

file(GLOB configs_src "configs/*.cpp")

file(GLOB rules_src "rules/*.cpp")

rdkit_library(CIPLabeler
              CIPMol.cpp CIPLabeler.cpp
              Mancude.cpp Digraph.cpp Node.cpp Edge.cpp Sort.cpp
              ${configs_src} ${rules_src}
              LINK_LIBRARIES GraphMol)

rdkit_headers(CIPLabeler.h
              DEST GraphMol/CIPLabeler)

rdkit_catch_test(testCIPLabeler catch_tests.cpp
                 CIPMol.cpp CIPLabeler.cpp
                 Mancude.cpp Digraph.cpp Node.cpp Edge.cpp Sort.cpp
                 ${configs_src} ${rules_src}
                 LINK_LIBRARIES SmilesParse GraphMol)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()









