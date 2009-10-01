add_library(ChemReactions Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp )

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES ChemReactions FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
