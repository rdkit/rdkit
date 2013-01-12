rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	      ReactionWriter.cpp ReactionDepict.cpp
              LINK_LIBRARIES Depictor FileParsers SubstructMatch ChemTransforms)

rdkit_headers(Reaction.h
              ReactionParser.h
              ReactionPickler.h DEST GraphMol/ChemReactions)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES
ChemReactions ChemTransforms Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
