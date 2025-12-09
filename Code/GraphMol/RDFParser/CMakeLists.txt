
rdkit_library(RDFParser
        RDFParser.cpp
        LINK_LIBRARIES GenericGroups Depictor GraphMol ChemReactions )
target_compile_definitions(RDFParser PRIVATE RDKIT_RDFPARSER_BUILD)

rdkit_headers(RDFParser.h.h
        DEST GraphMol/RDFParser)

rdkit_test(testRDFToMol testRDFToMol.cpp
        LINK_LIBRARIES MarvinParser FileParsers SubstructMatch )


