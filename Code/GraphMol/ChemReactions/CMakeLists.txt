add_library(ChemReactions Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp )

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES ChemReactions GraphMol FileParsers Depictor SmilesParse SubstructMatch RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
