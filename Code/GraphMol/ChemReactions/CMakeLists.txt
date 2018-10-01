if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    ADD_DEFINITIONS("-DRDK_USE_BOOST_SERIALIZATION")
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
      message("== Making EnumerateLibrary without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5)
  if(Boost_VERSION LESS  106200)
    set_property(SOURCE Enumerate/Enumerate.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -fpermissive ")
  endif()
endif()


remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_CHEMREACTIONS_BUILD)
rdkit_library(ChemReactions
              Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp
	            ReactionWriter.cpp ReactionDepict.cpp ReactionFingerprints.cpp ReactionUtils.cpp MoleculeParser.cpp ReactionRunner.cpp PreprocessRxn.cpp SanitizeRxn.cpp
                    Enumerate/Enumerate.cpp
                    Enumerate/EnumerationPickler.cpp
                    Enumerate/EvenSamplePairs.cpp

              LINK_LIBRARIES
              FilterCatalog Descriptors Fingerprints DataStructs Depictor
              FileParsers SubstructMatch ChemTransforms ${RDKit_SERIALIZATION_LIBS})

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
 ChemReactions ChemTransforms Descriptors Fingerprints
Subgraphs DataStructs  FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_SERIALIZATION_LIBS})

rdkit_test(testReactionFingerprints testReactionFingerprints.cpp LINK_LIBRARIES
ChemReactions Descriptors Fingerprints Subgraphs DataStructs ChemTransforms FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_SERIALIZATION_LIBS})

rdkit_test(testEnumeration Enumerate/testEnumerate.cpp LINK_LIBRARIES
ChemReactions ChemTransforms Descriptors Fingerprints Subgraphs DataStructs FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_SERIALIZATION_LIBS} )

rdkit_catch_test(rxnTestCatch catch_tests.cpp LINK_LIBRARIES ChemReactions ChemTransforms Descriptors Fingerprints
Subgraphs DataStructs  FileParsers SmilesParse SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_SERIALIZATION_LIBS} )

add_subdirectory(Wrap)
