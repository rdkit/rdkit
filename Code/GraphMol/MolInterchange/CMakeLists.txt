
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0")
    downloadAndCheckMD5("https://github.com/Tencent/rapidjson/archive/v1.1.0.tar.gz"
          "${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0.tar.gz"
          "badd12c511e081fec6c89c43a7027bce")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0.tar.gz
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/External)
    # comment out a line which causes a compilation error on some platforms
    # (based on the change which has already been applied to the RapidJSON master branch, see
    # https://github.com/Tencent/rapidjson/blob/ab1842a2dae061284c0a62dca1cc6d5e7e37e346/include/rapidjson/document.h#L414)
    file(READ ${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0/include/rapidjson/document.h RAPIDJSON_DOCUMENT_H)
    string(REGEX REPLACE
           "( *)(GenericStringRef& operator=\\(const GenericStringRef& rhs\\) { s = rhs\\.s. length = rhs\\.length. })" "\\1//\\2"
           RAPIDJSON_DOCUMENT_H "${RAPIDJSON_DOCUMENT_H}")
    file(WRITE ${CMAKE_SOURCE_DIR}/External/rapidjson-1.1.0/include/rapidjson/document.h "${RAPIDJSON_DOCUMENT_H}")
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
