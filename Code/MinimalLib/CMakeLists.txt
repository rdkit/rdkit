add_definitions(-DRDKIT_GRAPHMOL_BUILD)
rdkit_library(RDKit_Minimal
              minilib.cpp
              LINK_LIBRARIES SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral
                 )


rdkit_test(testMiniLib testMinilib.cpp LINK_LIBRARIES RDKit_Minimal_static)


