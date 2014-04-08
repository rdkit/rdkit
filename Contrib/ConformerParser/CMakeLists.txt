rdkit_library(ConformerParser ConformerParser.cpp
              LINK_LIBRARIES GraphMol)
              
rdkit_headers(ConformerParser.h DEST ../Contrib/ConformerParser)

rdkit_test(testConformerParser test.cpp 
           LINK_LIBRARIES ConformerParser
FileParsers SmilesParse
GraphMol RDGeometryLib RDGeneral)

add_subdirectory(Wrap)