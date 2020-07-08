remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLENUMERATOR_BUILD)
rdkit_library(MolEnumerator
              MolEnumerator.cpp PositionVariation.cpp LinkNode.cpp
              LINK_LIBRARIES FileParsers ChemReactions ChemTransforms 
              SubstructMatch GraphMol RDGeneral)

rdkit_headers(MolEnumerator.h DEST GraphMol/MolEnumerator)

rdkit_catch_test(testMolEnumerator catch_main.cpp enumerator_catch.cpp 
LINK_LIBRARIES MolEnumerator SmilesParse FileParsers )

#if(RDK_BUILD_PYTHON_WRAPPERS)
#add_subdirectory(Wrap)
#endif()
