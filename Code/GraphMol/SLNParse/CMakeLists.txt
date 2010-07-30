FIND_PACKAGE(BISON)
FIND_PACKAGE(FLEX)

find_package(Boost 1.39.0 COMPONENTS regex REQUIRED)
if(MSVC)
ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
# as of cmake 2.8.0 and boost 1.43 we need to do this to avoid
# forcing the user to build both static and dynamic versions
# of the regex library: 
ADD_DEFINITIONS("/D BOOST_ALL_DYN_LINK")  
endif()

FLEX_TARGET(SlnL sln.ll 
            ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysln.cpp 
            COMPILE_FLAGS "-Pyysln_")

BISON_TARGET(SlnY sln.yy 
             ${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.cpp 
             COMPILE_FLAGS "-pyysln_")

ADD_FLEX_BISON_DEPENDENCY(SlnL SlnY)

rdkit_library(SLNParse
              SLNParse.cpp SLNAttribs.cpp 
              ${BISON_SlnY_OUTPUTS} ${FLEX_SlnL_OUTPUTS}
              LINK_LIBRARIES GraphMol ${Boost_REGEX_LIBRARY} )

rdkit_test(testSLNParse test.cpp 
           LINK_LIBRARIES SLNParse SmilesParse SubstructMatch
GraphMol RDGeometryLib RDGeneral  ${Boost_REGEX_LIBRARY}
)

add_subdirectory(Wrap)


