find_package(Boost 1.39.0 COMPONENTS serialization REQUIRED)
message("==> ${Boost_SERIALIZATION_LIBRARY}")
rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	            ReactionWriter.cpp ReactionDepict.cpp ReactionFingerprints.cpp ReactionUtils.cpp MoleculeParser.cpp ReactionRunner.cpp PreprocessRxn.cpp SanitizeRxn.cpp
                    Enumerate/Enumerate.cpp
                    Enumerate/EnumerationPickler.cpp
                    Enumerate/EvenSamplePairs.cpp
                    
              LINK_LIBRARIES ${Boost_SERIALIZATION_LIBRARY}
                             FilterCatalog Descriptors Fingerprints DataStructs Depictor FileParsers SubstructMatch ChemTransforms)

rdkit_headers(Reaction.h
              ReactionParser.h
              ReactionPickler.h 
              ReactionFingerprints.h 
              ReactionUtils.h 
              ReactionRunner.h
              PreprocessRxn.h
              SanitizeRxn.h
              Enumerate/Enumerate.h
              Enumerate/EnumerateBase.h
              Enumerate/EnumerationPickler.h
              Enumerate/EnumerationStrategyBase.h
              Enumerate/CartesianProduct.h
              Enumerate/RandomSample.h
              Enumerate/RandomSampleAllBBs.h
              DEST GraphMol/ChemReactions)

rdkit_test(testReaction testReaction.cpp LINK_LIBRARIES
${Boost_SERIALIZATION_LIBRARY} ChemReactions ChemTransforms Descriptors Fingerprints Subgraphs DataStructs Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

rdkit_test(testReactionFingerprints testReactionFingerprints.cpp LINK_LIBRARIES
${Boost_SERIALIZATION_LIBRARY} ChemReactions Descriptors Fingerprints Subgraphs DataStructs ChemTransforms Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

rdkit_test(testEnumeration Enumerate/testEnumerate.cpp LINK_LIBRARIES
${Boost_SERIALIZATION_LIBRARY} ChemReactions ChemTransforms Descriptors Fingerprints Subgraphs DataStructs Depictor FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
