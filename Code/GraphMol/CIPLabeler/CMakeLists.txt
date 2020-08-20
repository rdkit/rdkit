
set(configs_src
    configs/Configuration.cpp
    configs/Sp2Bond.cpp
    configs/Tetrahedral.cpp)

set(rules_src
    rules/Rule1a.cpp
    rules/Rule1b.cpp
    rules/Rule2.cpp
    rules/Rule3.cpp
    rules/Rule4a.cpp
    rules/Rule4b.cpp
    rules/Rule4c.cpp
    rules/Rule5.cpp
    rules/Rule5New.cpp
    rules/Rule6.cpp
    rules/SequenceRule.cpp)

rdkit_library(CIPLabeler
              CIPMol.cpp CIPLabeler.cpp
              Mancude.cpp Digraph.cpp Node.cpp Edge.cpp Sort.cpp
              ${configs_src} ${rules_src}
              LINK_LIBRARIES GraphMol)
target_compile_definitions(CIPLabeler PRIVATE RDKIT_CIPLABELER_BUILD)

rdkit_headers(CIPLabeler.h
              DEST CIPLabeler)

rdkit_catch_test(testCIPLabeler catch_tests.cpp
                 CIPMol.cpp CIPLabeler.cpp
                 Mancude.cpp Digraph.cpp Node.cpp Edge.cpp Sort.cpp
                 ${configs_src} ${rules_src}
                 LINK_LIBRARIES SmilesParse GraphMol)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()









