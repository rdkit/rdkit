rdkit_library(ConformerParser ConformerParser
              LINK_LIBRARIES GraphMol)

rdkit_headers(ConformerParser.h DEST GraphMol/ConformerParser)

rdkit_test(testConformerParser test.cpp 
           LINK_LIBRARIES ConformerParser
FileParsers SmilesParse 
GraphMol RDGeometryLib RDGeneral 
)

add_subdirectory(Wrap)