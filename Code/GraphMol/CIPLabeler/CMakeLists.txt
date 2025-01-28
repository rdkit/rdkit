
set(configs_src
    configs/Configuration.cpp
    configs/Sp2Bond.cpp
    configs/Tetrahedral.cpp
    configs/AtropisomerBond.cpp)

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
              TooManyNodesException.h
              DEST GraphMol/CIPLabeler)

rdkit_catch_test(testCIPLabeler catch_tests.cpp 
                 CIPMol.cpp CIPLabeler.cpp
                 Mancude.cpp Digraph.cpp Node.cpp Edge.cpp Sort.cpp
                 ${configs_src} ${rules_src}
                 LINK_LIBRARIES SmilesParse GraphMol FileParsers MarvinParser)

# Set the BUILD macro to keep dll linkage consistency, since, in order
# to test private, unexported functions, we are building the test from
# sources, and NOT linking to the lib.
if(RDK_BUILD_CPP_TESTS)
    target_compile_definitions(testCIPLabeler PRIVATE RDKIT_CIPLABELER_BUILD)
endif()

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
