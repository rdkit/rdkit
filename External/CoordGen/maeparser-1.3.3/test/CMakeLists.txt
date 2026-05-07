include_directories(..)

find_package(Boost COMPONENTS filesystem iostreams unit_test_framework REQUIRED)

add_executable(unittest MainTestSuite.cpp BufferTest.cpp MaeBlockTest.cpp MaeParserTest.cpp
               ReaderTest.cpp WriterTest.cpp UsageDemo.cpp)

if(MAEPARSER_BUILD_SHARED_LIBS)
    target_compile_definitions(unittest PRIVATE "BOOST_ALL_DYN_LINK")
else(MAEPARSER_BUILD_SHARED_LIBS)
    target_compile_definitions(unittest PRIVATE "STATIC_MAEPARSER")
endif(MAEPARSER_BUILD_SHARED_LIBS)

get_filename_component(TEST_SAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
target_compile_definitions(unittest PRIVATE "TEST_SAMPLES_PATH=\"${TEST_SAMPLES_PATH}\"")

target_link_libraries(unittest maeparser ${Boost_LIBRARIES})

add_test(NAME unittest COMMAND ${CMAKE_CURRENT_BINARY_DIR}/unittest
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(MSVC AND MAEPARSER_BUILD_SHARED_LIBS)
    add_custom_command(TARGET unittest POST_BUILD COMMAND ${CMAKE_COMMAND}
    -E copy_if_different "${PROJECT_BINARY_DIR}/\$\(Configuration\)/maeparser.dll"
    "${CMAKE_CURRENT_BINARY_DIR}/\$\(Configuration\)/maeparser.dll")
endif(MSVC AND MAEPARSER_BUILD_SHARED_LIBS)

