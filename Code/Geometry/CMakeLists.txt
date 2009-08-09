add_library(RDGeometry point.cpp Transform2D.cpp Transform3D.cpp UniformGrid3D.cpp GridUtils.cpp )

add_executable(testTransforms testTransforms.cpp)
target_link_libraries(testTransforms  RDGeometry RDGeneral DataStructs)
add_test(testTransforms ${EXECUTABLE_OUTPUT_PATH}/testTransforms)

add_executable(testGrid testGrid.cpp)
target_link_libraries(testGrid  RDGeometry RDGeneral DataStructs)
add_test(testGrid ${EXECUTABLE_OUTPUT_PATH}/testGrid)

add_subdirectory(Wrap)
