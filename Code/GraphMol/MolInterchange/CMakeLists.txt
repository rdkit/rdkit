
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0")
    downloadAndCheckMD5("https://github.com/Tencent/rapidjson/archive/v1.1.0.tar.gz"
          "${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0.tar.gz"
          "badd12c511e081fec6c89c43a7027bce")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0.tar.gz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/External)
else()
  message("-- Found RapidJSON source in ${CMAKE_SOURCE_DIR}/External")
endif()

include_directories(${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0/include)


rdkit_library(MolInterchange
              Parser.cpp Writer.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(MolInterchange PRIVATE RDKIT_MOLINTERCHANGE_BUILD)

rdkit_headers(MolInterchange.h details.h
              DEST GraphMol/MolInterchange)

rdkit_test(molInterchangeTest1 test1.cpp
           LINK_LIBRARIES MolInterchange FileParsers SmilesParse)

rdkit_catch_test(molInterchangeCatchTest molinterchange_catch.cpp 
LINK_LIBRARIES MolInterchange FileParsers SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
