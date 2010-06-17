rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	      ReactionWriter.cpp ReactionDepict.cpp
              LINK_LIBRARIES Depictor FileParsers SubstructMatch)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES ChemReactions Depictor FileParsers SmilesParse SubstructMatch 
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
