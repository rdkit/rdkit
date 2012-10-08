FIND_PACKAGE(BISON)
FIND_PACKAGE(FLEX)

if(MSVC)
ADD_DEFINITIONS("/D YY_NO_UNISTD_H")
endif()

if(FLEX_EXECUTABLE)
  FLEX_TARGET(SmilesL smiles.ll 
              ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmiles.cpp 
             COMPILE_FLAGS "-Pyysmiles_" )
  FLEX_TARGET(SmartsL smarts.ll 
              ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmarts.cpp 
              COMPILE_FLAGS "-Pyysmarts_" )
endif(FLEX_EXECUTABLE)

if(BISON_EXECUTABLE)
  BISON_TARGET(SmilesY smiles.yy 
               ${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.cpp 
               COMPILE_FLAGS "-pyysmiles_" )
  BISON_TARGET(SmartsY smarts.yy 
               ${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.cpp 
               COMPILE_FLAGS "-pyysmarts_" )
  ADD_FLEX_BISON_DEPENDENCY(SmilesL SmilesY)
  ADD_FLEX_BISON_DEPENDENCY(SmartsL SmartsY)
endif(BISON_EXECUTABLE)

FILE(GLOB FLEX_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lex.*.cpp")
FILE(GLOB BISON_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tab.?pp")

get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
        "${extra_clean_files};${FLEX_OUTPUT_FILES};${BISON_OUTPUT_FILES}")


rdkit_library(SmilesParse
              SmilesParse.cpp SmilesParseOps.cpp 
              SmilesWrite.cpp SmartsWrite.cpp
              ${BISON_OUTPUT_FILES}
              ${FLEX_OUTPUT_FILES}
              LINK_LIBRARIES GraphMol)

rdkit_headers(primes.h
              SmartsWrite.h
              SmilesParse.h
              SmilesParseOps.h
              SmilesWrite.h DEST GraphMol/SmilesParse)

rdkit_test(smiTest1 test.cpp LINK_LIBRARIES SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(smaTest1 smatest.cpp LINK_LIBRARIES SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )

