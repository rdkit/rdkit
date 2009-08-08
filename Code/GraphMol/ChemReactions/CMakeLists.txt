add_library(ChemReactions Reaction.cpp MDLParser.cpp DaylightParser.cpp ReactionPickler.cpp )

add_executable(testReaction testReaction.cpp)
target_link_libraries(testReaction  ChemReactions FileParsers Depictor SmilesParse SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(testReaction ${EXECUTABLE_OUTPUT_PATH}/testReaction )

