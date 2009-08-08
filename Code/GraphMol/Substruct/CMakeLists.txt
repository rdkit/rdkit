add_library(SubstructMatch SubstructMatch.cpp SubstructUtils.cpp)

add_executable(substructMatchTest1 test1.cpp)
target_link_libraries(substructMatchTest1  SubstructMatch GraphMol RDGeneral RDGeometry )
add_test(substructMatchTest1 ${EXECUTABLE_OUTPUT_PATH}/substructMatchTest1 )

