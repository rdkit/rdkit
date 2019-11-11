remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_SCAFFOLDNETWORK_BUILD)
rdkit_library(ScaffoldNetwork
              ScaffoldNetwork.cpp
              LINK_LIBRARIES MolStandardize ChemReactions ChemTransforms SmilesParse SubstructMatch GraphMol RDGeneral)

rdkit_headers(ScaffoldNetwork.h DEST GraphMol/ScaffoldNetwor)

rdkit_catch_test(testScaffoldNetwork catch_tests.cpp 
LINK_LIBRARIES ScaffoldNetwork SmilesParse GraphMol RDGeneral 
${RDKit_THREAD_LIBS} )

#if(RDK_BUILD_PYTHON_WRAPPERS)
#add_subdirectory(Wrap)
#endif()
