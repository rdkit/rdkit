rdkit_python_extension(rdfragcatalog
                       rdfragcatalog.cpp FragCatalog.cpp FragCatParams.cpp 
                       FragCatGenerator.cpp FragFPGenerator.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
                       FragCatalog 
Subgraphs SubstructMatch SmilesParse Catalogs FileParsers GraphMol DataStructs RDGeometryLib RDGeneral 
 RDBoost) 

add_pytest(pyFragCatalog ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
