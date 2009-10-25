rdkit_python_extension(rdChemReactions
                       rdChemReactions.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
ChemReactions FileParsers SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost)

add_test(pyChemReactions python 
         ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py)

