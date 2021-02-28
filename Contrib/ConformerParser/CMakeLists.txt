rdkit_library(ConformerParser ConformerParser.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(SmilesParse PRIVATE RDKIT_CONFORMERPARSER_BUILD)

rdkit_test(testConformerParser test.cpp
           LINK_LIBRARIES ConformerParser FileParsers SmilesParse
           GraphMol RDGeometryLib RDGeneral)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()