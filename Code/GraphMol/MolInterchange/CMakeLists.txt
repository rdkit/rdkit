
if(NOT EXISTS "${RDKit_ExternalDir}/rapidjson/rapidjson-1.1.0")
    downloadAndCheckMD5("https://github.com/Tencent/rapidjson/archive/v1.1.0.tar.gz"
          "${RDKit_ExternalDir}/rapidjson/rapidjson-1.1.0.tar.gz"
          "badd12c511e081fec6c89c43a7027bce")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${RDKit_ExternalDir}/rapidjson/rapidjson-1.1.0.tar.gz
      WORKING_DIRECTORY ${RDKit_ExternalDir}/rapidjson)
    # backport some master branch improvements to avoid a compilation error on some platforms
    # (based on the change which has already been applied to the RapidJSON master branch, see
    # https://github.com/Tencent/rapidjson/blob/24b5e7a8b27f42fa16b96fc70aade9106cf7102f/include/rapidjson/document.h#L414
    # and https://patch-diff.githubusercontent.com/raw/Tencent/rapidjson/pull/1698.diff
    execute_process(COMMAND git --git-dir= apply --ignore-whitespace
      ${RDKit_ExternalDir}/rapidjson/patch-rapidjson-1.1.0.diff
      WORKING_DIRECTORY ${RDKit_ExternalDir}/rapidjson)
else()
    message("-- Found RapidJSON source in ${RDKit_ExternalDir}/rapidjson")
endif()

include_directories(${RDKit_ExternalDir}/rapidjson/rapidjson-1.1.0/include)


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
