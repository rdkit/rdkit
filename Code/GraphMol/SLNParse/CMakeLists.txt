if(RDK_USE_FLEXBISON)
  FIND_PACKAGE(BISON)
  FIND_PACKAGE(FLEX)
else(RDK_USE_FLEXBISON)
  set(FLEX_EXECUTABLE "")
  set(BISON_EXECUTABLE "")
endif(RDK_USE_FLEXBISON)

set(T_LIBS ${Boost_LIBRARIES})
find_package(Boost 1.39.0 COMPONENTS regex REQUIRED)
set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})

if(MSVC)
  ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
  # as of cmake 2.8.0 and boost 1.43 we need to do this to avoid
  # forcing the user to build both static and dynamic versions
  # of the regex library:
  ADD_DEFINITIONS("/D BOOST_ALL_DYN_LINK") 
endif(MSVC)

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
              LINK_LIBRARIES GraphMol ${Boost_REGEX_LIBRARY} )

rdkit_headers(SLNAttribs.h
              SLNParse.h
              SLNParseOps.h DEST GraphMol/SLNParse)

rdkit_test(testSLNParse test.cpp 
           LINK_LIBRARIES SLNParse SmilesParse SubstructMatch
GraphMol RDGeometryLib RDGeneral  ${Boost_REGEX_LIBRARY}
)

add_subdirectory(Wrap)


