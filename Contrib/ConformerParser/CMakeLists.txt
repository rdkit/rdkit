rdkit_library(ConformerParser ConformerParser.cpp
              LINK_LIBRARIES GraphMol)

rdkit_test(testConformerParser test.cpp 
           LINK_LIBRARIES ConformerParser
FileParsers SmilesParse
GraphMol RDGeometryLib RDGeneral)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()