
rdkit_library(MarvinParser
              MarvinParser.cpp
              MarvinWriter.cpp
              MarvinDefs.cpp
              LINK_LIBRARIES GenericGroups Depictor GraphMol ChemReactions )
target_compile_definitions(MarvinParser PRIVATE RDKIT_MARVINPARSER_BUILD)

rdkit_headers(MarvinParser.h
              DEST GraphMol/MarvinParse)

rdkit_test(testMrvToMol testMrvToMol.cpp
           LINK_LIBRARIES MarvinParser FileParsers SubstructMatch )


