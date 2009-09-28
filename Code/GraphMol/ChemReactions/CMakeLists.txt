rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
              LINK_LIBRARIES FileParsers SubstructMatch)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES ChemReactions)

add_subdirectory(Wrap)
