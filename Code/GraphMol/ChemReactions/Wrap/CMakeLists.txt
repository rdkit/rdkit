rdkit_python_extension(rdChemReactions
  rdChemReactions.cpp
  DEST Chem 
  LINK_LIBRARIES ChemReactions SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometry DataStructs RDGeneral RDBoost )
add_test( pyChemReactions python ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py )

