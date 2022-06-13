rdkit_library(ScaffoldNetwork
              ScaffoldNetwork.cpp
              LINK_LIBRARIES MolStandardize ChemReactions ChemTransforms SmilesParse 
              SubstructMatch GraphMol RDGeneral)
target_compile_definitions(ScaffoldNetwork PRIVATE RDKIT_SCAFFOLDNETWORK_BUILD)

rdkit_headers(ScaffoldNetwork.h DEST GraphMol/ScaffoldNetwork)

rdkit_catch_test(testScaffoldNetwork catch_tests.cpp 
LINK_LIBRARIES ScaffoldNetwork SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
