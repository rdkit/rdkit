rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	            ReactionWriter.cpp ReactionDepict.cpp ReactionFingerprints.cpp ReactionUtils.cpp MoleculeParser.cpp ReactionRunner.cpp
              LINK_LIBRARIES Descriptors Fingerprints DataStructs Depictor FileParsers SubstructMatch ChemTransforms)

rdkit_headers(Reaction.h
              ReactionParser.h
              ReactionPickler.h 
              ReactionFingerprints.h 
              ReactionUtils.h 
              ReactionRunner.h DEST GraphMol/ChemReactions)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES
ChemReactions ChemTransforms Descriptors Fingerprints Subgraphs DataStructs Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

rdkit_test(testReactionFingerprints testReactionFingerprints.cpp LINK_LIBRARIES
ChemReactions Descriptors Fingerprints Subgraphs DataStructs ChemTransforms Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
