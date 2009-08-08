add_library(Depictor RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp)

add_executable(testDepictor testDepictor.cpp)
target_link_libraries(testDepictor Depictor SmilesParse FileParsers GraphMol RDGeneral RDGeometry )
add_test(testDepictor ${EXECUTABLE_OUTPUT_PATH}/testDepictor )


