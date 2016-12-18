rdkit_python_extension(rdChemReactions
                       Enumerate.cpp
                       rdChemReactions.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
ChemReactions FilterCatalog ChemTransforms Descriptors Fingerprints Subgraphs DataStructs Depictor FileParsers SmilesParse SubstructMatch GraphMol Catalogs FilterCatalog RDGeneral RDGeometryLib RDBoost ${Boost_SERIALIZATION_LIBRARY})

add_pytest(pyChemReactions
         ${CMAKE_CURRENT_SOURCE_DIR}/testReactionWrapper.py)

add_pytest(pyChemReactionEnumerations
         ${CMAKE_CURRENT_SOURCE_DIR}/testEnumerations.py)

add_pytest(pyChemReactionSanitize
         ${CMAKE_CURRENT_SOURCE_DIR}/testSanitize.py)
