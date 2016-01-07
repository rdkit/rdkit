rdkit_python_extension(rdChemReactions
                       Enumerate.cpp
                       rdChemReactions.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
ChemReactions ChemTransforms Descriptors Fingerprints Subgraphs DataStructs Depictor FileParsers SmilesParse SubstructMatch GraphMol Catalogs FilterCatalog RDGeneral RDGeometryLib RDBoost)

add_pytest(pyChemReactions
         ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py)

add_pytest(pyChemReactionEnumerations
         ${CMAKE_CURRENT_SOURCE_DIR}/testEnumerations.py)
