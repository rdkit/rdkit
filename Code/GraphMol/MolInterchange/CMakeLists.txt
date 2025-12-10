
rdkit_library(MolInterchange
              Parser.cpp Writer.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(MolInterchange PRIVATE RDKIT_MOLINTERCHANGE_BUILD)
target_compile_definitions(MolInterchange PRIVATE BOOST_JSON_STATIC_LINK)

rdkit_headers(MolInterchange.h details.h
              DEST GraphMol/MolInterchange)

rdkit_test(molInterchangeTest1 test1.cpp
           LINK_LIBRARIES MolInterchange FileParsers SmilesParse)

rdkit_catch_test(molInterchangeCatchTest molinterchange_catch.cpp 
LINK_LIBRARIES MolInterchange FileParsers SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
