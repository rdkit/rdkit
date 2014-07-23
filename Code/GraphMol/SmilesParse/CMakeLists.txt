if(RDK_USE_FLEXBISON)
  FIND_PACKAGE(BISON)
  FIND_PACKAGE(FLEX)
else(RDK_USE_FLEXBISON)
  set(FLEX_EXECUTABLE "")
  set(BISON_EXECUTABLE "")
endif(RDK_USE_FLEXBISON)

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
  SET(FLEX_OUTPUT_FILES ${FLEX_SmilesL_OUTPUTS} ${FLEX_SmartsL_OUTPUTS})
else(FLEX_EXECUTABLE)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmiles.cpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmiles.cpp COPYONLY)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmarts.cpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/lex.yysmarts.cpp COPYONLY)
  FILE(GLOB FLEX_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lex.*.cpp")	
endif(FLEX_EXECUTABLE)

if(BISON_EXECUTABLE)
  BISON_TARGET(SmilesY smiles.yy 
               ${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.cpp 
               COMPILE_FLAGS "-pyysmiles_" )
  BISON_TARGET(SmartsY smarts.yy 
               ${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.cpp 
               COMPILE_FLAGS "-pyysmarts_" )
  SET(BISON_OUTPUT_FILES ${BISON_SmilesY_OUTPUTS} ${BISON_SmartsY_OUTPUTS})
  if(FLEX_EXECUTABLE)
    ADD_FLEX_BISON_DEPENDENCY(SmilesL SmilesY)
    ADD_FLEX_BISON_DEPENDENCY(SmartsL SmartsY)
  endif(FLEX_EXECUTABLE)
else(BISON_EXECUTABLE)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.hpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.hpp COPYONLY)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.cpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/smiles.tab.cpp COPYONLY)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.hpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.hpp COPYONLY)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.cpp.cmake 
                 ${CMAKE_CURRENT_SOURCE_DIR}/smarts.tab.cpp COPYONLY)
  FILE(GLOB BISON_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tab.?pp")
endif(BISON_EXECUTABLE)

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

rdkit_test(smiTest1 test.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib )

rdkit_test(smaTest1 smatest.cpp LINK_LIBRARIES SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib )

