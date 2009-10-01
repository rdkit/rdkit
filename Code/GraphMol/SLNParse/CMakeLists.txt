FIND_PACKAGE(BISON)
FIND_PACKAGE(FLEX)

if(MSVC)
ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
endif()

FLEX_TARGET(SlnL sln.ll 
            ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysln.cpp 
            COMPILE_FLAGS "-P yysln_")

BISON_TARGET(SlnY sln.yy 
             ${CMAKE_CURRENT_SOURCE_DIR}/sln.tab.cpp 
             COMPILE_FLAGS "-p yysln_")

ADD_FLEX_BISON_DEPENDENCY(SlnL SlnY)

rdkit_library(SLNParse
              InputFiller.cpp SLNParse.cpp SLNAttribs.cpp 
              ${BISON_SlnY_OUTPUTS} ${FLEX_SlnL_OUTPUTS}
              LINK_LIBRARIES GraphMol)

rdkit_test(testSLNParse test.cpp 
           LINK_LIBRARIES SLNParse SmilesParse SubstructMatch
GraphMol RDGeometryLib RDGeneral
)

add_subdirectory(Wrap)


