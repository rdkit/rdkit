rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	      ReactionWriter.cpp ReactionDepict.cpp
              LINK_LIBRARIES Depictor FileParsers SubstructMatch)

rdkit_headers(Reaction.h
              ReactionParser.h
              ReactionPickler.h DEST GraphMol/ChemReactions)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES ChemReactions_static Depictor_static FileParsers_static SmilesParse_static SubstructMatch_static
GraphMol_static RDGeneral_static RDGeometryLib_static )

add_subdirectory(Wrap)
