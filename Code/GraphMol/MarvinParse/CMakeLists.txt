
rdkit_library(MarvinParser
              MarvinParser.cpp
              MarvinWriter.cpp
              MarvinDefs.cpp
              ${maesupplier}
              LINK_LIBRARIES GenericGroups Depictor GraphMol ChemReactions ${RDK_MAEPARSER_LIBS} ${regex_lib} ${link_iostreams} ${zlib_lib})
target_compile_definitions(MarvinParser PRIVATE RDKIT_MARVINPARSER_BUILD)

rdkit_headers(MarvinParser.h
              DEST GraphMol/MarvinParser)

rdkit_test(testMrvToMol testMrvToMol.cpp
           LINK_LIBRARIES MarvinParser FileParsers SubstructMatch )


