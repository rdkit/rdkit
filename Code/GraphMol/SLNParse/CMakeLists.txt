if(RDK_USE_FLEXBISON)
  FIND_PACKAGE(BISON)
  FIND_PACKAGE(FLEX)
else(RDK_USE_FLEXBISON)
  set(FLEX_EXECUTABLE "")
  set(BISON_EXECUTABLE "")
endif(RDK_USE_FLEXBISON)

if(MSVC)
  ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
endif(MSVC)

if(CMAKE_COMPILER_IS_CLANG AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 15)
# the output from bison triggers some warnings we can safely ignore with clang 16+
set_source_files_properties(sln.tab.cpp PROPERTIES COMPILE_OPTIONS "-Wno-unused-but-set-variable")
endif()



if(FLEX_EXECUTABLE)
  FLEX_TARGET(SlnL sln.ll
              ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysln.cpp
              COMPILE_FLAGS "-Pyysln_")
  SET(FLEX_OUTPUT_FILES ${FLEX_SlnL_OUTPUTS})
else(FLEX_EXECUTABLE)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/lex.yysln.cpp.cmake
                 ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysln.cpp COPYONLY)
  FILE(GLOB FLEX_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lex.*.cpp")
endif(FLEX_EXECUTABLE)

if(BISON_EXECUTABLE)
  BISON_TARGET(SlnY sln.yy
               ${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.cpp
               COMPILE_FLAGS "-pyysln_")
  SET(BISON_OUTPUT_FILES ${BISON_SlnY_OUTPUTS})
  if(FLEX_EXECUTABLE)
    ADD_FLEX_BISON_DEPENDENCY(SlnL SlnY)
  endif(FLEX_EXECUTABLE)
else(BISON_EXECUTABLE)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.hpp.cmake
                 ${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.hpp COPYONLY)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.cpp.cmake
                 ${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.cpp COPYONLY)
  FILE(GLOB BISON_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tab.?pp")
endif(BISON_EXECUTABLE)



rdkit_library(SLNParse
              SLNParse.cpp SLNAttribs.cpp
              ${BISON_OUTPUT_FILES} ${FLEX_OUTPUT_FILES}
              LINK_LIBRARIES GraphMol RDGeneral )
target_compile_definitions(SLNParse PRIVATE RDKIT_SLNPARSE_BUILD)

rdkit_headers(SLNAttribs.h
              SLNParse.h
              SLNParseOps.h DEST GraphMol/SLNParse)

rdkit_test(testSLNParse test.cpp
           LINK_LIBRARIES SLNParse SmilesParse SubstructMatch )

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()
