rdkit_library(Abbreviations
              Abbreviations.cpp
              AbbreviationsUtils.cpp
              LINK_LIBRARIES SmilesParse 
              SubstructMatch GraphMol RDGeneral)
target_compile_definitions(Abbreviations PRIVATE RDKIT_ABBREVIATIONS_BUILD)

rdkit_headers(Abbreviations.h DEST GraphMol/Abbreviations)

rdkit_catch_test(testAbbreviations catch_tests.cpp 
LINK_LIBRARIES Abbreviations SmilesParse FileParsers )

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()
