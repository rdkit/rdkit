rdkit_python_extension(rdChemReactions
                       rdChemReactions.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
ChemReactions Depictor FileParsers SmilesParse SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost)

add_pytest(pyChemReactions
         ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py)

