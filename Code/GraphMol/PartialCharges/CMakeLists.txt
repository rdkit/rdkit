
rdkit_library(PartialCharges
              GasteigerCharges.cpp GasteigerParams.cpp
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(PartialCharges PRIVATE RDKIT_PARTIALCHARGES_BUILD)

rdkit_headers(GasteigerCharges.h
              GasteigerParams.h DEST GraphMol/PartialCharges)

rdkit_test(testPartialCharges test1.cpp LINK_LIBRARIES  PartialCharges
SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
