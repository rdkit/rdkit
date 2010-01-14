FIND_PACKAGE(BISON)
FIND_PACKAGE(FLEX)

if(MSVC)
ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
endif()


FLEX_TARGET(SmilesL smiles.ll 
            ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmiles.cpp 
            COMPILE_FLAGS "-Pyysmiles_" )
FLEX_TARGET(SmartsL smarts.ll 
            ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmarts.cpp 
            COMPILE_FLAGS "-Pyysmarts_" )
BISON_TARGET(SmilesY smiles.yy 
             ${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.cpp 
             COMPILE_FLAGS "-pyysmiles_" )
BISON_TARGET(SmartsY smarts.yy 
             ${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.cpp 
             COMPILE_FLAGS "-pyysmarts_" )

ADD_FLEX_BISON_DEPENDENCY(SmilesL SmilesY)
ADD_FLEX_BISON_DEPENDENCY(SmartsL SmartsY)

rdkit_library(SmilesParse
              InputFiller.cpp SmilesParse.cpp SmilesParseOps.cpp 
              SmilesWrite.cpp SmartsWrite.cpp
              ${BISON_SmilesY_OUTPUTS} ${BISON_SmartsY_OUTPUTS} 
              ${FLEX_SmilesL_OUTPUTS}  ${FLEX_SmartsL_OUTPUTS}
              LINK_LIBRARIES GraphMol)

rdkit_test(smiTest1 test.cpp LINK_LIBRARIES SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(smaTest1 smatest.cpp LINK_LIBRARIES SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )

